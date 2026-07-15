#!/usr/bin/env python3
"""
Transition-robustness screen for the airfoil shortlists.

Re-evaluates each element's top candidates under three boundary-layer
transition environments:

  clean  n_crit=9              XFoil/NeuralFoil default: clean wind tunnel.
                               Matches the main screen - the baseline.
  track  n_crit=6              realistic outdoor/gusty ambient turbulence.
  abuse  xtr=0.05 (forced)     leading edge tripped: rubber pickup, bugs,
                               rain. No laminar flow survives.

A section whose CLmax / L/D collapses in the track or abuse case was
depending on laminar flow it won't get on lap 15 of endurance.

NOTE ON SCOPE: airfoil_screen.py now GATES on the abuse case (a section must
reach cl_target + margin when tripped, or it is disqualified outright). So the
shortlist this script reads has already survived that filter. What this script
adds is the *intermediate* `track` case and the full clean->track->abuse slope,
which show HOW a section degrades rather than merely whether it clears the bar.
Run it to choose between survivors, not to find failures - the screen already
removed those.

Evaluated at the LOW-Re end of each element's band by default (the worst
corner: laminar separation bubbles are nastiest at low Re). --re design
reproduces the old design-point behavior.

Usage (after a screen run):
  python robustness.py                     # reads screen_results/
  python robustness.py --results DIR --top 8 --re design
Outputs:
  <results>/robustness_<profile>.png       slope charts across the 3 cases
  <results>/robustness_summary.csv         all numbers + retention ratios

CAUTION ON COLUMN NAMES: CLmax_clean here is measured at ONE Re (the low-Re
end, or design Re with --re design). The screen's `<profile>_CLmax` column in
results.xlsx is the WORST CASE ACROSS the whole Re band. They are different
quantities and will not match. Do not diff them.
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import (
    ALPHA, CASES, CONFIG, PROFILES, ld_at, load_kulfan, load_run_config,
    masked_polar,
)

HIGHLIGHT = ["#d62728", "#1f77b4", "#2ca02c"]


def metrics(aero, cl_target):
    """CLmax and L/D at target CL from one polar.

    Uses common.masked_polar, so it applies the SAME confidence threshold and
    the SAME sweep-boundary guard as the main screen. (This script used to use
    a 0.5 valid-fraction where the screen used 0.6, and took argmax with no
    boundary check - so a polar whose CLmax sat on the edge of the sweep, i.e.
    was never actually resolved, would be reported as a real measurement.)

    L/D is np.nan - NOT 0.0 - when the section cannot reach cl_target. Those
    are different failures and the caller must not confuse them.
    """
    m = masked_polar(aero, ALPHA)
    if m is None:
        return np.nan, np.nan
    a, cl, cd, cm, i_max = m
    return float(cl[i_max]), ld_at(cl, cd, i_max, cl_target)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="screen_results")
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--re", choices=["lo", "design"], default="lo",
                    help="Re at which to run the sweep (default: lo = worst corner)")
    args = ap.parse_args()

    import neuralfoil as nf

    # Analyze at the fidelity the results were PRODUCED at, not whatever the
    # defaults happen to be today.
    load_run_config(args.results)
    kulfan = load_kulfan(args.results)
    print(f"model_size={CONFIG['model_size']} (from run_config.json)  "
          f"Re={args.re}  cases={list(CASES)}")

    xlsx = os.path.join(args.results, "results.xlsx")
    plotdir = os.path.join(args.results, "plots")   # all figures live in plots/
    os.makedirs(plotdir, exist_ok=True)
    rows = []

    for pname, prof in PROFILES.items():
        try:
            top = pd.read_excel(xlsx, sheet_name=f"{pname}_top").head(args.top)
        except ValueError:
            continue
        re_eval = prof["re_list"][0 if args.re == "lo" else 1]
        # Same per-element TE blunt the screen applied. Without it we would be
        # re-analyzing the as-drawn knife edge, not the section that was ranked.
        te_needed = CONFIG["te_min_mm"] / (prof["chord"] * 1000.0)
        data = {}
        for _, r in top.iterrows():
            kp = kulfan.get(r["name"])
            if kp is None:
                print(f"  ! {r['name']} missing from kulfan_params.json - skipped")
                continue
            kp = dict(kp)
            kp["TE_thickness"] = max(float(kp["TE_thickness"]), te_needed)
            per_case = {}
            for cname, kw in CASES.items():
                aero = nf.get_aero_from_kulfan_parameters(
                    kulfan_parameters=kp, alpha=ALPHA, Re=re_eval,
                    model_size=CONFIG["model_size"], **kw)
                per_case[cname] = metrics(aero, prof["cl_target"])
            data[r["name"]] = per_case
            clean, track, abuse = (per_case[c] for c in CASES)
            rows.append(dict(
                profile=pname, name=r["name"], Re=re_eval,
                CLmax_clean=clean[0], CLmax_track=track[0], CLmax_abuse=abuse[0],
                LD_clean=clean[1], LD_track=track[1], LD_abuse=abuse[1],
                # Explicit, rather than encoding "unreachable" as a 0.0 that
                # looks like a measurement and poisons any later mean().
                target_unreachable_abuse=bool(np.isnan(abuse[1])),
                CLmax_retention=(abuse[0] / clean[0]
                                 if clean[0] and not np.isnan(clean[0]) else np.nan),
                LD_retention=(abuse[1] / clean[1]
                              if clean[1] and not np.isnan(clean[1]) else np.nan),
            ))

        if not data:
            print(f"  ! {pname}: no candidates to plot")
            continue

        # slope charts: one line per airfoil across the 3 cases
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
        xs = np.arange(len(CASES))
        for k, (name, pc) in enumerate(data.items()):
            style = (dict(color=HIGHLIGHT[k], lw=2, marker="o", zorder=3)
                     if k < 3 else dict(color="0.65", lw=1, marker="o",
                                        alpha=0.6, zorder=1))
            lbl = f"{k+1}. {name}" if k < 3 else None
            ld = [pc[c][1] for c in CASES]
            ax1.plot(xs, [pc[c][0] for c in CASES], label=lbl, **style)
            ax2.plot(xs, ld, label=lbl, **style)
            # NaN = target unreachable. matplotlib just drops the point, which
            # would silently hide the single most important failure mode, so
            # mark it.
            for xi, v in zip(xs, ld):
                if np.isnan(v):
                    ax2.plot(xi, 0, marker="x", ms=9, mew=2,
                             color=style.get("color", "0.65"), zorder=4)
        for ax, ttl in ((ax1, "CLmax"), (ax2, f"L/D at CL {prof['cl_target']}")):
            ax.set_xticks(xs)
            ax.set_xticklabels([f"{c}\n{'ncrit 9' if c=='clean' else 'ncrit 6' if c=='track' else 'tripped LE'}"
                                for c in CASES], fontsize=9)
            ax.set_title(ttl, fontsize=11)
            ax.grid(alpha=0.3)
        ax2.text(0.02, 0.02, "x = cannot reach target CL", transform=ax2.transAxes,
                 fontsize=8, color="0.35")
        ax1.legend(fontsize=8)
        fig.suptitle(f"{pname.upper()} transition robustness  (Re {re_eval//1000}k)",
                     fontsize=12)
        fig.tight_layout()
        out = os.path.join(plotdir, f"robustness_{pname}.png")
        fig.savefig(out, dpi=130)
        plt.close(fig)
        print(f"wrote {out}")

    df = pd.DataFrame(rows)
    csv = os.path.join(args.results, "robustness_summary.csv")
    df.to_csv(csv, index=False)
    print(f"wrote {csv}")

    if df.empty:
        return
    print("\nMost robust per profile (highest L/D retention clean -> abuse):")
    for pname in df.profile.unique():
        sub = df[df.profile == pname].sort_values("LD_retention", ascending=False)
        sub = sub[~sub.LD_retention.isna()]
        if sub.empty:
            print(f"  {pname:14s} - no candidate reaches target CL when tripped")
            continue
        r = sub.iloc[0]
        print(f"  {pname:14s} {r['name']:<20} keeps {r.LD_retention*100:4.0f}% L/D, "
              f"{r.CLmax_retention*100:4.0f}% CLmax when tripped")

    dead = df[df.target_unreachable_abuse]
    if len(dead):
        print(f"\n! {len(dead)} shortlisted section(s) cannot reach target CL with a "
              f"FULLY tripped LE: {', '.join(dead.name.unique())}")
        print("  This is NOT a disqualification. The screen gates on the racing")
        print("  condition (n_crit 6, light LE grit); the fully-tripped case is rain /")
        print("  heavy contamination. These sections race fine and are fragile in the")
        print("  wet - decide whether that trade is acceptable for this element.")


if __name__ == "__main__":
    main()
