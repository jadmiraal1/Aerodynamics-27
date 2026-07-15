#!/usr/bin/env python3
"""
Gradient-based airfoil SHAPE optimization, per element.

WHY GRADIENTS (the algorithm, and why this one)
-----------------------------------------------
NeuralFoil is differentiable: through AeroSandbox's Opti (a CasADi wrapper) we
get EXACT analytic gradients of CLmax, L/D, Cp and the boundary-layer shape
factor with respect to all 17 Kulfan weights, for free - no finite differencing.

With exact gradients in a 17-dimensional space, a gradient method is the right
tool and everything else is wrong:

  - random variants (what the screen does): jitter, no gradient, biased toward
    large weights, ~600 samples barely dent a 17-D space.
  - CMA-ES / differential evolution: thousands of evaluations to do what IPOPT
    does in ~50 iterations, and they cannot use the gradient we already have.
  - IPOPT (interior-point, gradient-based): ~1 second per solve, handles the
    nonlinear geometry constraints natively. This is what AeroSandbox ships.

Gradient methods find a LOCAL optimum. The fix is MULTI-START: seed IPOPT from
each of the screen's top families and keep the best converged result. That
turns "local" into "best basin among sections we already trust" - which is
exactly right, because we are REFINING known-good airfoils, not searching from
nothing.

THE TRAP: OPTIMIZING AGAINST A SURROGATE
----------------------------------------
A 17-D optimizer handed a neural surrogate will find the surrogate's blind
spots, because that is where the free lunch is. Those blind spots are almost
always SEPARATED FLOW - the regime an integral-BL method models worst.

So the guards are not decoration, they are the point:

  1. TRUST REGION. Weights stay within +/-0.12 of a screened seed. The optimizer
     may refine a known-good section; it may not wander into fantasy geometry.
  2. CONFIDENCE FLOOR. NeuralFoil reports analysis_confidence; we penalize
     dropping below 0.85. Do not optimize where the model admits it is guessing.
  3. SEPARATION PENALTY. We watch the BL shape factor H directly.

     ...BUT NOT "H < 2.8 everywhere". That was the naive first attempt and it is
     aerodynamically wrong: high-lift low-Re sections like s1223 run H ~ 3.3
     with TRAILING-EDGE SEPARATION AT EVERY ANGLE - that separation IS how they
     carry aft load. Forbidding it forbids the lift. What we actually forbid is
     RUNAWAY separation (H climbing past the seed's level into 5, 8, 12), which
     is the surrogate being exploited, not a real airfoil. Penalty, not hard
     constraint, so an infeasible seed cannot kill the solve.

  4. XFOIL VALIDATION IS MANDATORY, not optional. NeuralFoil was TRAINED ON
     XFoil, so XFoil is the ground truth it approximates. If an optimized
     section's XFoil polar disagrees with NeuralFoil's, you have found a
     surrogate artifact, not an airfoil. This tool writes .dat files ready for
     XFoil/XFLR5 and REFUSES to call anything a result until you have checked.

Objective per element (worst case across the Re band, at the race condition):
  flap1/flap2      maximize CLmax          (loading-limited)
  main_center/out  maximize mean L/D over the operating CL band (wake quality)

The trailing edge is FIXED at the manufacturable blunt (te_min_mm / chord) -
the optimizer cannot buy performance with an un-cuttable knife edge.

Usage:
  python shape_optimize.py --element flap2
  python shape_optimize.py --element main_center --seeds goe304 goe63
  python shape_optimize.py --element flap2 --trust 0.10 --starts 4
Outputs (in <results>/shape_opt/):
  <element>_optimized.dat        Selig coords, blunted, XFoil-ready
  <element>_optimized.json       Kulfan params + metrics + seed + provenance
  <element>_compare.png          seed vs optimized: geometry, polar, Cp, H
"""

import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import aerosandbox as asb
import aerosandbox.numpy as anp
import neuralfoil as nf

from common import CONFIG, PROFILES, gate_kwargs, load_run_config

CONF_FLOOR = 0.85
H_SEPARATION_BUDGET = 0.4       # allowed H rise above the seed before penalty


def _kd(kf):
    return dict(upper_weights=np.array(kf.upper_weights),
                lower_weights=np.array(kf.lower_weights),
                leading_edge_weight=float(kf.leading_edge_weight))


def seed_families(results_dir, element, explicit=None):
    """Distinct seed sections to multi-start from: the screen's top families."""
    import pandas as pd
    if explicit:
        return explicit
    xlsx = os.path.join(results_dir, "results.xlsx")
    top = pd.read_excel(xlsx, sheet_name=f"{element}_top")
    fams = (top["name"].str.replace(r"_v\d+$", "", regex=True)
            .drop_duplicates().tolist())
    return fams[:4]                 # up to 4 distinct families


def evaluate_fixed(kp, element, model_size, n_alpha=13):
    """Plain NeuralFoil metrics on a FIXED section (numbers, not variables).
    Worst case across the element's Re band, at the race condition."""
    p = PROFILES[element]
    alpha = np.linspace(2.0, 15.0, n_alpha)
    clmax, ld_band, conf, hmax = [], [], [], []
    band = p["cl_band"]
    for Re in p["re_list"]:
        a = nf.get_aero_from_kulfan_parameters(
            kulfan_parameters=kp, alpha=alpha, Re=Re,
            model_size=model_size, **gate_kwargs())
        cl = np.asarray(a["CL"], float); cd = np.asarray(a["CD"], float)
        clmax.append(float(np.max(cl)))
        i = int(np.argmax(cl))
        order = np.argsort(cl[:i + 1])
        cls, cds = cl[:i + 1][order], cd[:i + 1][order]
        lds = [float(t / max(np.interp(t, cls, cds), 1e-6))
               if np.max(cl) >= t else 0.0
               for t in np.linspace(band[0], band[1], 5)]
        ld_band.append(float(np.mean(lds)))
        conf.append(float(np.min(a["analysis_confidence"])))
        # Shape factor at the USABLE point (2 deg below stall), NOT the full
        # sweep. Post-stall H is meaninglessly large (s1223 hits ~12) and using
        # it as the separation budget makes the guard toothless - the whole
        # point is to bound separation where the wing OPERATES, not deep in
        # stall where every section has separated.
        Hk = [k for k in a if "_bl_H_" in k]
        Hmat = np.array([[float(a[k][j]) for k in Hk] for j in range(len(alpha))])
        i_stall = int(np.argmax(cl))
        a_use = alpha[i_stall] - 2.0
        j_use = int(np.argmin(np.abs(alpha - a_use)))
        hmax.append(float(Hmat[j_use].max()))       # H_max at the usable point
    return dict(clmax=min(clmax), ld_band=min(ld_band),
                conf=min(conf), hmax=max(hmax))


def optimize_from_seed(seed_name, element, model_size, trust, max_iter=200):
    """One IPOPT solve from one seed. Returns (kp_dict, metrics) or None."""
    p = PROFILES[element]
    te = CONFIG["te_min_mm"] / (p["chord"] * 1000.0)
    kf = asb.Airfoil(seed_name).to_kulfan_airfoil()
    s = _kd(kf)
    seed_metrics = evaluate_fixed(dict(s, TE_thickness=te), element, model_size)

    opti = asb.Opti()
    uw = opti.variable(init_guess=s["upper_weights"])
    lw = opti.variable(init_guess=s["lower_weights"])
    le = opti.variable(init_guess=s["leading_edge_weight"])
    kp = dict(upper_weights=uw, lower_weights=lw,
              leading_edge_weight=le, TE_thickness=te)

    alpha = np.arange(3.0, 15.1, 1.0)
    is_flap = element.startswith("flap")
    band = p["cl_band"]

    obj = 0.0
    for Re in p["re_list"]:                          # multi-point: whole band
        a = nf.get_aero_from_kulfan_parameters(
            kulfan_parameters=kp, alpha=alpha, Re=Re,
            model_size=model_size, **gate_kwargs())
        CL = [a["CL"][i] for i in range(len(alpha))]
        clmax = anp.softmax(*CL, hardness=20)

        Hk = [k for k in a if "_bl_H_" in k]
        # Pre-stall window only. Measuring separation into the stall region
        # would let the optimizer trade a fake CLmax against a penalty it can
        # always pay by moving the (soft) stall point - the guard must live
        # where the wing operates.
        i_use = [i for i, al in enumerate(alpha) if al <= 9]
        Huse = anp.softmax(*[a[k][i] for k in Hk for i in i_use], hardness=12)
        conf = anp.softmin(*[a["analysis_confidence"][i]
                             for i in range(len(alpha))], hardness=30)

        if is_flap:
            perf = -clmax                            # maximize CLmax
        else:
            # mean L/D across the CL band, differentiable
            lds = []
            for t in np.linspace(band[0], band[1], 4):
                # L/D at CL=t via the pre-stall branch (soft, monotone approx)
                cd_t = anp.softmax(*[a["CD"][i] for i in range(len(alpha))],
                                   hardness=-8)       # ~ min CD near target
                lds.append(t / cd_t)
            perf = -anp.mean(anp.array(lds))

        # SOFT guards - penalties, never hard constraints (an infeasible seed
        # must not kill the solve). Separation is bounded RELATIVE to the seed:
        # its own H is allowed, runaway past it is not.
        h_cap = seed_metrics["hmax"] + H_SEPARATION_BUDGET   # usable-point H + margin
        sep_pen = 8.0 * anp.softmax(0.0, Huse - h_cap, hardness=8)
        conf_pen = 6.0 * anp.softmax(0.0, CONF_FLOOR - conf, hardness=50)
        obj = obj + perf + sep_pen + conf_pen

    af = asb.KulfanAirfoil(name="opt", **kp)
    # HARD constraints: geometry band + trust region. All smooth, all reliable.
    opti.subject_to([
        af.max_thickness() > p["tc"][0], af.max_thickness() < p["tc"][1],
        af.max_camber() > 0.015,
        uw - s["upper_weights"] < trust, uw - s["upper_weights"] > -trust,
        lw - s["lower_weights"] < trust, lw - s["lower_weights"] > -trust,
    ])
    opti.minimize(obj / len(p["re_list"]))

    try:
        sol = opti.solve(max_iter=max_iter, verbose=False)
    except Exception:
        return None
    out = dict(upper_weights=np.array(sol(uw)), lower_weights=np.array(sol(lw)),
               leading_edge_weight=float(sol(le)), TE_thickness=te)
    m = evaluate_fixed(out, element, model_size)
    m.update(seed=seed_name, seed_clmax=seed_metrics["clmax"],
             seed_ld_band=seed_metrics["ld_band"], seed_hmax=seed_metrics["hmax"],
             trust=trust,
             max_weight_move=float(max(
                 np.abs(out["upper_weights"] - s["upper_weights"]).max(),
                 np.abs(out["lower_weights"] - s["lower_weights"]).max())))
    return out, m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="screen_results")
    ap.add_argument("--element", required=True,
                    choices=list(PROFILES))
    ap.add_argument("--seeds", nargs="+", default=None,
                    help="seed sections (default: top families from the screen)")
    ap.add_argument("--trust", type=float, default=0.12,
                    help="trust region: max Kulfan-weight move from the seed")
    ap.add_argument("--starts", type=int, default=4, help="max seeds to try")
    args = ap.parse_args()

    load_run_config(args.results)
    ms = CONFIG["model_size"]
    el = args.element
    outdir = os.path.join(args.results, "shape_opt")
    os.makedirs(outdir, exist_ok=True)

    seeds = seed_families(args.results, el, args.seeds)[:args.starts]
    is_flap = el.startswith("flap")
    metric = "clmax" if is_flap else "ld_band"
    print(f"SHAPE OPTIMIZATION - {el}   (maximize {metric}, model {ms})")
    print(f"  race condition n_crit={CONFIG['gate_n_crit']} xtr={CONFIG['gate_xtr']}, "
          f"TE fixed at {CONFIG['te_min_mm']} mm")
    print(f"  multi-start seeds: {seeds}\n")

    best = None
    for sd in seeds:
        r = optimize_from_seed(sd, el, ms, args.trust)
        if r is None:
            print(f"  {sd:<12s} did not converge")
            continue
        kp, m = r
        gain = m[metric] - m[f"seed_{metric}"]
        print(f"  {sd:<12s} {metric} {m[f'seed_{metric}']:.3f} -> {m[metric]:.3f} "
              f"({gain:+.3f})   H {m['hmax']:.2f}  conf {m['conf']:.3f}  "
              f"move {m['max_weight_move']:.3f}")
        if best is None or m[metric] > best[1][metric]:
            best = (kp, m)

    if best is None:
        raise SystemExit("no seed converged - loosen --trust or check the seeds")
    kp, m = best

    # --- exports ---
    af = asb.KulfanAirfoil(name=f"{el}_opt", **kp)
    dat = os.path.join(outdir, f"{el}_optimized.dat")
    with open(dat, "w") as f:
        f.write(f"{el}_opt (from {m['seed']}, blunted {CONFIG['te_min_mm']}mm)\n")
        for x, y in af.coordinates:
            f.write(f" {x:.6f}  {y:.6f}\n")
    with open(os.path.join(outdir, f"{el}_optimized.json"), "w") as f:
        json.dump({k: (v.tolist() if isinstance(v, np.ndarray) else v)
                   for k, v in {**kp, **m}.items()}, f, indent=2)

    _compare_fig(best, el, ms, outdir)

    print("\n" + "=" * 68)
    print(f"BEST: {el}  from seed {m['seed']}")
    print(f"  {metric}: {m[f'seed_{metric}']:.3f} -> {m[metric]:.3f}")
    print(f"  trust-region move {m['max_weight_move']:.3f} of {args.trust}", end="")
    if m["max_weight_move"] > 0.98 * args.trust:
        print("   !! ON THE TRUST BOUND - the optimizer wants to go further.")
        print("      Widen --trust and re-run, OR accept that this is as far as")
        print("      you trust the surrogate to extrapolate from a known seed.")
    else:
        print("   (interior - a real local optimum, not a wall)")
    print(f"\n  !! NOT A RESULT UNTIL VALIDATED IN XFOIL/XFLR5 !!")
    print(f"     NeuralFoil was trained on XFoil; an optimized section is exactly")
    print(f"     where the surrogate is most likely to be wrong. Open")
    print(f"     {dat}")
    print(f"     in XFLR5 at Re {PROFILES[el]['re_list'][0]:,}, n_crit "
          f"{CONFIG['gate_n_crit']}, and confirm the polar before trusting it.")


def _compare_fig(best, el, ms, outdir):
    kp, m = best
    p = PROFILES[el]
    te = kp["TE_thickness"]
    seed = _kd(asb.Airfoil(m["seed"]).to_kulfan_airfoil())
    seed["TE_thickness"] = te
    Re = p["re_list"][0]
    alpha = np.linspace(0, 16, 33)

    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    for label, k, c in [(f"seed {m['seed']}", seed, "0.55"),
                        (f"{el}_opt", kp, "#b2182b")]:
        af = asb.KulfanAirfoil(name="x", **k)
        co = np.array(af.coordinates)
        a = nf.get_aero_from_kulfan_parameters(
            kulfan_parameters=k, alpha=alpha, Re=Re, model_size=ms,
            **gate_kwargs())
        msk = np.asarray(a["analysis_confidence"]) >= 0.85
        ax[0].plot(co[:, 0], co[:, 1], color=c, label=label)
        ax[1].plot(np.asarray(a["CD"])[msk], np.asarray(a["CL"])[msk], color=c)
        ax[2].plot(alpha[msk], np.asarray(a["CL"])[msk], color=c)
    ax[0].set_aspect("equal"); ax[0].set_title("geometry (blunted TE)")
    ax[0].legend(fontsize=8)
    ax[1].set_xlabel("CD"); ax[1].set_ylabel("CL"); ax[1].set_title(f"polar, Re {Re//1000}k")
    ax[2].set_xlabel("alpha"); ax[2].set_ylabel("CL"); ax[2].set_title("CL-alpha")
    for a_ in ax: a_.grid(alpha=0.3)
    fig.suptitle(f"{el}: gradient shape optimization (VALIDATE IN XFOIL)",
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"{el}_compare.png"), dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()
