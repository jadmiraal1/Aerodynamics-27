#!/usr/bin/env python3
"""
Transition-robustness screen for the airfoil shortlists.

Re-evaluates each element's top candidates under three boundary-layer
transition environments at the element's design-point Re:

  clean  n_crit=9              XFoil/NeuralFoil default: clean wind tunnel.
                               Matches the main screen - the baseline.
  track  n_crit=6              realistic outdoor/gusty ambient turbulence.
  abuse  xtr=0.05 (forced)     leading edge tripped: rubber pickup, bugs,
                               rain. No laminar flow survives.

A section whose CLmax / L/D collapses in the track or abuse case was
depending on laminar flow it won't get on lap 15 of endurance.

Usage (after a screen run):
  python robustness.py                     # reads screen_results/
  python robustness.py --results DIR --top 8
Outputs:
  <results>/robustness_<profile>.png       slope charts across the 3 cases
  <results>/robustness_summary.csv         all numbers + retention ratios
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from airfoil_screen import PROFILES, ALPHA, MODEL_SIZE, CONF_MIN

CASES = {
    "clean": dict(n_crit=9.0),
    "track": dict(n_crit=6.0),
    "abuse": dict(n_crit=9.0, xtr_upper=0.05, xtr_lower=0.05),
}
HIGHLIGHT = ["#d62728", "#1f77b4", "#2ca02c"]


def metrics(aero, cl_target):
    """CLmax and L/D at target CL from one polar (confidence-masked)."""
    conf = np.asarray(aero["analysis_confidence"], float)
    m = conf >= CONF_MIN
    if m.mean() < 0.5:
        return np.nan, np.nan
    cl = np.asarray(aero["CL"], float)[m]
    cd = np.asarray(aero["CD"], float)[m]
    i = int(np.argmax(cl))
    clmax = float(cl[i])
    pre_cl, pre_cd = cl[: i + 1], cd[: i + 1]
    o = np.argsort(pre_cl)
    if clmax < cl_target:
        return clmax, 0.0
    cd_t = float(np.interp(cl_target, pre_cl[o], pre_cd[o]))
    return clmax, cl_target / max(cd_t, 1e-6)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="screen_results")
    ap.add_argument("--top", type=int, default=8)
    args = ap.parse_args()

    import aerosandbox as asb
    import neuralfoil as nf

    xlsx = os.path.join(args.results, "results.xlsx")
    datd = os.path.join(args.results, "shortlist_dat")
    rows = []

    for pname, prof in PROFILES.items():
        try:
            top = pd.read_excel(xlsx, sheet_name=f"{pname}_top").head(args.top)
        except ValueError:
            continue
        re_design = prof["re_list"][1]
        data = {}
        for _, r in top.iterrows():
            dat = os.path.join(datd, f"{r['name']}.dat")
            if not os.path.exists(dat):
                continue
            af = asb.Airfoil(name=r["name"], coordinates=np.loadtxt(dat, skiprows=1))
            per_case = {}
            for cname, kw in CASES.items():
                aero = nf.get_aero_from_airfoil(
                    af, alpha=ALPHA, Re=re_design, model_size=MODEL_SIZE, **kw)
                per_case[cname] = metrics(aero, prof["cl_target"])
            data[r["name"]] = per_case
            clean, track, abuse = (per_case[c] for c in CASES)
            rows.append(dict(
                profile=pname, name=r["name"],
                CLmax_clean=clean[0], CLmax_track=track[0], CLmax_abuse=abuse[0],
                LD_clean=clean[1], LD_track=track[1], LD_abuse=abuse[1],
                CLmax_retention=abuse[0] / clean[0] if clean[0] else np.nan,
                LD_retention=abuse[1] / clean[1] if clean[1] else np.nan,
            ))

        # slope charts: one line per airfoil across the 3 cases
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
        xs = np.arange(len(CASES))
        for k, (name, pc) in enumerate(data.items()):
            style = (dict(color=HIGHLIGHT[k], lw=2, marker="o", zorder=3)
                     if k < 3 else dict(color="0.65", lw=1, marker="o",
                                        alpha=0.6, zorder=1))
            lbl = f"{k+1}. {name}" if k < 3 else None
            ax1.plot(xs, [pc[c][0] for c in CASES], label=lbl, **style)
            ax2.plot(xs, [pc[c][1] for c in CASES], label=lbl, **style)
        for ax, ttl in ((ax1, "CLmax"), (ax2, f"L/D at CL {prof['cl_target']}")):
            ax.set_xticks(xs)
            ax.set_xticklabels([f"{c}\n{'ncrit 9' if c=='clean' else 'ncrit 6' if c=='track' else 'tripped LE'}"
                                for c in CASES], fontsize=9)
            ax.set_title(ttl, fontsize=11)
            ax.grid(alpha=0.3)
        ax1.legend(fontsize=8)
        fig.suptitle(f"{pname.upper()} transition robustness  (Re {re_design//1000}k)",
                     fontsize=12)
        fig.tight_layout()
        out = os.path.join(args.results, f"robustness_{pname}.png")
        fig.savefig(out, dpi=130)
        plt.close(fig)
        print(f"wrote {out}")

    df = pd.DataFrame(rows)
    csv = os.path.join(args.results, "robustness_summary.csv")
    df.to_csv(csv, index=False)
    print(f"wrote {csv}")

    print("\nMost robust per profile (highest L/D retention clean -> abuse):")
    for pname in df.profile.unique():
        sub = df[df.profile == pname].sort_values("LD_retention", ascending=False)
        r = sub.iloc[0]
        print(f"  {pname:14s} {r['name']:<20} keeps {r.LD_retention*100:4.0f}% L/D, "
              f"{r.CLmax_retention*100:4.0f}% CLmax when tripped")


if __name__ == "__main__":
    main()
