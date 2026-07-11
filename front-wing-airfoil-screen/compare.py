#!/usr/bin/env python3
"""
Side-by-side comparison sheets for airfoil_screen.py results.

For each element profile, produces ONE figure overlaying the top candidates:
  - airfoil geometry (true aspect ratio)
  - CL vs alpha at the profile's design Re
  - drag polar (CL vs CD)
  - total score bar chart
Top 3 get color + labels; the rest are ghosted gray for context.

Usage (after a screen run):
  python compare.py                        # reads screen_results/
  python compare.py --results other_dir --top 10
Outputs: <results>/compare_<profile>.png
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from airfoil_screen import PROFILES, MODEL_SIZE, CONF_MIN

# wider sweep than the screen: shows the full drag bucket incl. the
# lower branch (these sections have zero-lift angles near -8 deg)
PLOT_ALPHA = np.arange(-10.0, 18.1, 0.5)

HIGHLIGHT = ["#d62728", "#1f77b4", "#2ca02c"]      # ranks 1-3
GHOST = dict(color="0.65", lw=0.9, alpha=0.55, zorder=1)


def load_dat(path):
    pts = np.loadtxt(path, skiprows=1)
    return pts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="screen_results")
    ap.add_argument("--top", type=int, default=8, help="candidates per sheet")
    args = ap.parse_args()

    import aerosandbox as asb
    import neuralfoil as nf

    xlsx = os.path.join(args.results, "results.xlsx")
    datd = os.path.join(args.results, "shortlist_dat")

    for pname, prof in PROFILES.items():
        try:
            top = pd.read_excel(xlsx, sheet_name=f"{pname}_top").head(args.top)
        except ValueError:
            continue
        re_design = prof["re_list"][1]          # middle = design-point Re

        fig, ax = plt.subplots(2, 2, figsize=(13, 8.5))
        (ax_geo, ax_cla), (ax_polar, ax_score) = ax

        names, scores = [], []
        for rank, (_, row) in enumerate(top.iterrows()):
            name = row["name"]
            dat = os.path.join(datd, f"{name}.dat")
            if not os.path.exists(dat):
                continue
            coords = load_dat(dat)
            af = asb.Airfoil(name=name, coordinates=coords)
            aero = nf.get_aero_from_airfoil(af, alpha=PLOT_ALPHA, Re=re_design,
                                            model_size=MODEL_SIZE)
            m = np.asarray(aero["analysis_confidence"]) >= CONF_MIN
            cl = np.asarray(aero["CL"])[m]
            cd = np.asarray(aero["CD"])[m]
            al = PLOT_ALPHA[m]

            if rank < 3:
                style = dict(color=HIGHLIGHT[rank], lw=2.0, zorder=3)
                geo_style = dict(color=HIGHLIGHT[rank], lw=1.6, zorder=3)
            else:
                style = GHOST
                geo_style = GHOST
            label = f"{rank+1}. {name}" if rank < 3 else None

            ax_geo.plot(coords[:, 0], coords[:, 1], label=label, **geo_style)
            ax_cla.plot(al, cl, label=label, **style)
            ax_polar.plot(cd, cl, label=label, **style)
            names.append(name)
            scores.append(row[f"score_{pname}"])

        # score bars (all candidates, rank order top->bottom)
        y = np.arange(len(names))
        colors = [HIGHLIGHT[i] if i < 3 else "0.65" for i in range(len(names))]
        ax_score.barh(y, scores, color=colors, alpha=0.85)
        ax_score.set_yticks(y)
        ax_score.set_yticklabels(names, fontsize=8)
        ax_score.invert_yaxis()
        ax_score.set_xlabel(f"score_{pname}")
        ax_score.set_xlim(min(scores) * 0.95, max(scores) * 1.02)

        # targets / reference lines
        ax_cla.axhline(prof["cl_target"], color="0.3", ls="--", lw=1)
        ax_cla.text(PLOT_ALPHA[0] + 0.3, prof["cl_target"] + 0.02, "CL target",
                    fontsize=8, color="0.3")
        ax_polar.axhline(prof["cl_target"], color="0.3", ls="--", lw=1)

        ax_geo.set_aspect("equal")
        ax_geo.set_title("geometry", fontsize=10)
        ax_cla.set_title(f"CL vs alpha  (Re {re_design//1000}k)", fontsize=10)
        ax_cla.set_xlabel("alpha [deg]"); ax_cla.set_ylabel("CL")
        ax_polar.set_title("drag polar", fontsize=10)
        ax_polar.set_xlabel("CD"); ax_polar.set_ylabel("CL")
        ax_score.set_title("total score", fontsize=10)
        for a_ in (ax_geo, ax_cla, ax_polar):
            a_.grid(alpha=0.3)
            a_.legend(fontsize=8, loc="best")

        fig.suptitle(
            f"{pname.upper()}  -  c={prof['chord']*1000:.0f}mm, "
            f"CL target {prof['cl_target']}, top {len(names)} of shortlist",
            fontsize=12)
        fig.tight_layout()
        out = os.path.join(args.results, f"compare_{pname}.png")
        fig.savefig(out, dpi=130)
        plt.close(fig)
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
