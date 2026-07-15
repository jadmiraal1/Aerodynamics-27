#!/usr/bin/env python3
"""
XFoil validation: does the NeuralFoil ranking survive the real solver?

WHY THIS IS NOT OPTIONAL
------------------------
NeuralFoil is a neural network TRAINED ON XFoil. So XFoil is the ground truth
it approximates, and every number in this whole pipeline - the screen, the
robustness sweep, the shape optimizer - is a surrogate estimate of what XFoil
would say. That is fine for ranking 2,000 airfoils fast. It is NOT fine as a
final answer, for two reasons:

  1. The surrogate has error, and the error is largest exactly where we care -
     near CLmax, at low Re, in separated flow.
  2. The shape optimizer actively SEEKS the surrogate's blind spots (that is
     where the free lunch is). An optimized section is the single most likely
     place for NeuralFoil and XFoil to disagree.

This stage runs XFoil at the SAME race condition the screen gated on (n_crit,
xtr, Re band) and plots NeuralFoil vs XFoil side by side. If they agree, you
have earned confidence. If they diverge, believe XFoil - and treat the
NeuralFoil ranking of that section as suspect.

REQUIRES THE XFOIL BINARY on PATH (or pass --xfoil-command path\\to\\xfoil.exe).
XFoil is free: https://web.mit.edu/drela/Public/web/xfoil/ . On Windows put
xfoil.exe next to this script or on PATH. Without it, use --mock to test the
harness (a second NeuralFoil model stands in - proves the plumbing, NOT the
physics).

Validates the BLUNTED geometry - the section that was screened and will be cut,
not the as-drawn one. (See DESIGN_JUSTIFICATION.md on the TE blunt.)

Usage:
  python validate.py --element flap2                 # top 3 of the shortlist
  python validate.py --element flap2 --optimized     # the shape-optimized one
  python validate.py --element flap2 --names s1223_v19 be6699_v04
  python validate.py --element flap2 --mock          # no xfoil, test the harness
Outputs (in <results>/validation/):
  <element>_validation.png     NeuralFoil vs XFoil, per section
  validation_summary.csv       agreement metrics, all sections
"""

import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import aerosandbox as asb
import neuralfoil as nf

from common import CONFIG, PROFILES, gate_kwargs, load_kulfan, load_run_config

ALPHA = np.arange(-2.0, 16.1, 0.5)
C_NF = "#1f77b4"       # NeuralFoil
C_XF = "#b2182b"       # XFoil (ground truth)


def blunt_kp(kp, element):
    """kp with the manufacturable TE, at this element's chord."""
    te = CONFIG["te_min_mm"] / (PROFILES[element]["chord"] * 1000.0)
    out = dict(kp)
    out["TE_thickness"] = max(float(out.get("TE_thickness", 0.0)), te)
    return out


def run_neuralfoil(kp, Re, model_size):
    a = nf.get_aero_from_kulfan_parameters(
        kulfan_parameters=kp, alpha=ALPHA, Re=Re, model_size=model_size,
        **gate_kwargs())
    m = np.asarray(a["analysis_confidence"]) >= CONFIG["conf_min"]
    return dict(alpha=ALPHA[m], CL=np.asarray(a["CL"])[m],
                CD=np.asarray(a["CD"])[m], CM=np.asarray(a["CM"])[m])


def run_xfoil(kp, Re, xfoil_command, mock_model=None):
    """XFoil at the race condition. mock_model != None -> NeuralFoil stand-in."""
    if mock_model is not None:
        # HARNESS TEST ONLY: a different NeuralFoil model pretending to be XFoil.
        # Proves the comparison/plot logic. Says NOTHING about real agreement.
        return run_neuralfoil(kp, Re, mock_model)

    af = asb.KulfanAirfoil(name="v", **kp)
    gk = gate_kwargs()
    xf = asb.XFoil(airfoil=af, Re=Re, n_crit=gk["n_crit"],
                   xtr_upper=gk["xtr_upper"], xtr_lower=gk["xtr_lower"],
                   xfoil_command=xfoil_command, max_iter=100, timeout=60)
    out = xf.alpha(ALPHA)                      # XFoil drops points it can't solve
    if "CL" not in out or len(out.get("alpha", [])) == 0:
        return None
    return dict(alpha=np.asarray(out["alpha"]), CL=np.asarray(out["CL"]),
                CD=np.asarray(out["CD"]),
                CM=np.asarray(out.get("CM", np.full_like(out["CL"], np.nan))))


def metrics(pol, cl_target):
    """CLmax and L/D at target from a polar dict."""
    if pol is None or len(pol["CL"]) < 3:
        return dict(clmax=np.nan, ld_target=np.nan, a_clmax=np.nan)
    cl, cd, al = pol["CL"], pol["CD"], pol["alpha"]
    i = int(np.argmax(cl))
    clmax = float(cl[i])
    ld = np.nan
    if clmax >= cl_target:
        o = np.argsort(cl[:i + 1])
        cd_t = float(np.interp(cl_target, cl[:i + 1][o], cd[:i + 1][o]))
        ld = cl_target / max(cd_t, 1e-6)
    return dict(clmax=clmax, ld_target=ld, a_clmax=float(al[i]))


def validate_section(name, kp, element, model_size, xfoil_command, mock):
    p = PROFILES[element]
    Re = p["re_list"][0]                       # low-Re end (worst corner)
    kpb = blunt_kp(kp, element)
    nf_pol = run_neuralfoil(kpb, Re, model_size)
    xf_pol = run_xfoil(kpb, Re, xfoil_command, mock_model=mock)

    nm = metrics(nf_pol, p["cl_target"])
    xm = metrics(xf_pol, p["cl_target"])
    d_clmax = nm["clmax"] - xm["clmax"]
    # relative L/D disagreement, the number that actually matters for ranking
    ld_rel = (abs(nm["ld_target"] - xm["ld_target"]) / xm["ld_target"]
              if np.isfinite(xm["ld_target"]) and xm["ld_target"] else np.nan)
    verdict = ("no XFoil" if xf_pol is None else
               "AGREE" if abs(d_clmax) < 0.05 and (np.isnan(ld_rel) or ld_rel < 0.10)
               else "DISAGREE")
    return dict(name=name, Re=Re,
                nf_clmax=nm["clmax"], xf_clmax=xm["clmax"], d_clmax=d_clmax,
                nf_ld=nm["ld_target"], xf_ld=xm["ld_target"], ld_rel_err=ld_rel,
                verdict=verdict), nf_pol, xf_pol


def figure(rows, pols, element, outdir, mock):
    n = len(rows)
    fig, ax = plt.subplots(n, 3, figsize=(13, 3.2 * n), squeeze=False)
    xf_label = "NeuralFoil (mock)" if mock else "XFoil"
    p = PROFILES[element]
    for r, (row, (nfp, xfp)) in enumerate(zip(rows, pols)):
        a_cl, a_po, a_ld = ax[r]
        for pol, c, lab in [(nfp, C_NF, "NeuralFoil"), (xfp, C_XF, xf_label)]:
            if pol is None or len(pol["CL"]) < 2:
                continue
            a_cl.plot(pol["alpha"], pol["CL"], color=c, label=lab)
            a_po.plot(pol["CD"], pol["CL"], color=c, label=lab)
            good = pol["CD"] > 1e-6
            a_ld.plot(pol["alpha"][good], pol["CL"][good] / pol["CD"][good], color=c)
        for a_ in (a_cl, a_po, a_ld):
            a_.grid(alpha=0.3)
        a_cl.axhline(p["cl_target"], color="0.4", ls="--", lw=0.8)
        a_po.axhline(p["cl_target"], color="0.4", ls="--", lw=0.8)
        a_cl.set_ylabel("CL"); a_cl.set_xlabel("alpha")
        a_po.set_xlabel("CD"); a_ld.set_xlabel("alpha"); a_ld.set_ylabel("L/D")
        a_cl.legend(fontsize=8, loc="lower right")
        tag = (f"{row['name']}   NF CLmax {row['nf_clmax']:.2f} vs "
               f"{'mock' if mock else 'XF'} {row['xf_clmax']:.2f}   "
               f"[{row['verdict']}]")
        a_cl.set_title(tag, fontsize=9, loc="left",
                       color=("#b2182b" if row["verdict"] == "DISAGREE" else "0.2"))
    banner = ("HARNESS TEST - 'XFoil' is a second NeuralFoil model, NOT real "
              "XFoil" if mock else
              f"NeuralFoil vs XFoil at the race condition "
              f"(Re {rows[0]['Re']//1000}k, n_crit {CONFIG['gate_n_crit']}, "
              f"xtr {CONFIG['gate_xtr']})")
    fig.suptitle(f"{element} validation — {banner}", fontweight="bold",
                 color=("#b2182b" if mock else "black"))
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    out = os.path.join(outdir, f"{element}_validation.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="screen_results")
    ap.add_argument("--element", required=True, choices=list(PROFILES) + ["all"])
    ap.add_argument("--names", nargs="+", default=None,
                    help="sections to validate (default: top 3 of the shortlist)")
    ap.add_argument("--optimized", action="store_true",
                    help="validate <element>_optimized from shape_optimize.py")
    ap.add_argument("--top", type=int, default=3)
    ap.add_argument("--xfoil-command", default="xfoil",
                    help="xfoil binary (e.g. xfoil.exe or a full path)")
    ap.add_argument("--mock", action="store_true",
                    help="no XFoil: a 2nd NeuralFoil model stands in (harness test)")
    ap.add_argument("--skip-if-no-xfoil", action="store_true",
                    help="exit cleanly (not an error) if XFoil is missing - for "
                         "the pipeline, so a run without XFoil does not fail")
    args = ap.parse_args()

    load_run_config(args.results)
    ms = CONFIG["model_size"]
    outdir = os.path.join(args.results, "validation")
    os.makedirs(outdir, exist_ok=True)
    mock = "small" if args.mock else None     # small != the screen model, on purpose

    if not mock and __import__("shutil").which(args.xfoil_command) is None:
        msg = (f"XFoil binary '{args.xfoil_command}' not found on PATH.\n"
               f"  Install it (free): https://web.mit.edu/drela/Public/web/xfoil/\n"
               f"  or pass --xfoil-command C:\\path\\to\\xfoil.exe\n"
               f"  or use --mock to test the harness without XFoil.")
        if args.skip_if_no_xfoil:
            print("VALIDATION skipped - " + msg.split("\n")[0])
            print("  (install XFoil to validate the NeuralFoil ranking; the "
                  "shortlist is a surrogate estimate until you do)")
            return
        raise SystemExit(msg)

    elements = list(PROFILES) if args.element == "all" else [args.element]
    for el in elements:
        _validate_one(el, args, ms, mock, outdir)


def _validate_one(el, args, ms, mock, outdir):

    # --- gather (name, kulfan) pairs to validate ---
    items = []
    if args.optimized:
        j = os.path.join(args.results, "shape_opt", f"{el}_optimized.json")
        if not os.path.exists(j):
            raise SystemExit(f"{j} not found - run shape_optimize.py first")
        d = json.load(open(j))
        kp = {k: (np.array(d[k]) if isinstance(d[k], list) else d[k])
              for k in ("upper_weights", "lower_weights",
                        "leading_edge_weight", "TE_thickness")}
        items.append((f"{el}_optimized", kp))
    else:
        kulfan = load_kulfan(args.results)
        names = args.names
        if names is None:
            top = pd.read_excel(os.path.join(args.results, "results.xlsx"),
                                sheet_name=f"{el}_top")
            names = top["name"].head(args.top).tolist()
        for nm in names:
            if nm in kulfan:
                items.append((nm, kulfan[nm]))
            else:
                print(f"  ! {nm} not in kulfan_params.json - skipped")

    print(f"VALIDATION - {el}   NeuralFoil({ms}) vs "
          f"{'MOCK' if mock else 'XFoil'}")
    print(f"  race condition: Re {PROFILES[el]['re_list'][0]:,}, "
          f"n_crit {CONFIG['gate_n_crit']}, xtr {CONFIG['gate_xtr']}\n")

    rows, pols = [], []
    for nm, kp in items:
        row, nfp, xfp = validate_section(nm, kp, el, ms, args.xfoil_command, mock)
        rows.append(row); pols.append((nfp, xfp))
        print(f"  {nm:<16s} CLmax NF {row['nf_clmax']:.3f} "
              f"{'MK' if mock else 'XF'} {row['xf_clmax']:.3f} "
              f"(d {row['d_clmax']:+.3f})   L/D err "
              f"{row['ld_rel_err']*100 if np.isfinite(row['ld_rel_err']) else float('nan'):4.0f}%"
              f"   -> {row['verdict']}")

    out = figure(rows, pols, el, outdir, mock)
    pd.DataFrame(rows).to_csv(
        os.path.join(outdir, "validation_summary.csv"), index=False)

    dis = [r for r in rows if r["verdict"] == "DISAGREE"]
    print(f"\n  wrote {out}")
    if mock:
        print("  (MOCK run - proves the harness only. Install XFoil for a real check.)")
    elif dis:
        print(f"  !! {len(dis)} section(s) DISAGREE with XFoil: "
              f"{', '.join(r['name'] for r in dis)}")
        print("     Believe XFoil. The NeuralFoil ranking of these is suspect -")
        print("     re-check before committing them to CFD or a mould.")
    else:
        print("  All validated sections agree with XFoil within tolerance.")


if __name__ == "__main__":
    main()
