#!/usr/bin/env python3
"""
Sensitivity of the airfoil selection to the scoring weights.

THE QUESTION
------------
The scoring weights in PROFILES are engineering JUDGMENT: "wake quality matters
0.40 for the center element, stall gentleness 0.25". Nobody derived those
numbers and nobody can. That is only a problem if the ANSWER depends on them.

If the same section wins across a wide spread of plausible weightings, the
judgment was not load-bearing and the selection is robust. If the winner flips
whenever a weight is nudged, then the weights - not the aerodynamics - are
choosing the wing. This script measures which of those two worlds we are in.

METHOD
------
1. MONTE CARLO. Weight vectors are drawn from a Dirichlet distribution centred
   on the nominal weights, w ~ Dir(w_nom * kappa). The concentration kappa sets
   the spread: high = draws stay near nominal, low = wide disagreement about
   what matters. Each draw re-scores the whole qualifier pool; the winner is
   recorded. The result is a WIN PROBABILITY per section.

   THE HEADLINE NUMBER IS THE FAMILY WIN PROBABILITY, NOT THE VARIANT.
   Rivals are almost always CST siblings of the winner (goe63_v16 vs goe63_v28
   vs goe63_v03). A mould is cut per FAMILY, not per variant. If the weights
   only reorder siblings, they are not choosing the wing - the aerodynamic
   decision was already made. Reporting the variant number alone understates
   robustness badly: flap1's winner takes 56% of draws at variant level, but
   s1223 takes 100% of them at family level.

2. ONE-AT-A-TIME SWEEP. Each weight is swept 0 -> 1 with the remainder
   renormalized, tracking the rank of the nominal winner. Identifies which
   single criterion, if any, the result actually hinges on.

3. SIMPLEX VERTICES. Score with each metric alone at weight 1.0. These are the
   extreme corners of the weight space. A section winning several corners cannot
   be dislodged by any interior weighting.

4. MARGIN. Score gap between rank 1 and rank 2. A winner with a negligible
   margin is not a winner regardless of the weights.

5. RANK STABILITY. Spearman correlation of each perturbed ranking against the
   nominal one, over the top 200 candidates. Tests whether the ORDER survives,
   not merely the winner.

Usage:
  python weight_sensitivity.py                       # reads screen_results/
  python weight_sensitivity.py --results DIR --draws 5000 --concentration 20

Outputs:
  <results>/plots/weight_sensitivity_<element>.png   one figure per element
  <results>/weight_sensitivity.csv                   win probabilities, all elements
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import PROFILES, load_run_config, normalized_metrics

METRICS = ["CLmax", "CL_usable", "stall_gentle", "Cm_low", "LD_max", "LD_band"]

# Muted, print-safe palette. No decorative colour.
C_WINNER = "#b2182b"
C_RIVAL = "#4393c3"
C_GHOST = "#cccccc"
C_GRID = "0.85"
C_TEXT = "0.30"

plt.rcParams.update({
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "axes.edgecolor": "0.4",
    "axes.linewidth": 0.8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "legend.frameon": False,
    "figure.dpi": 130,
})


def build_matrix(df, pname):
    """(names, M) where M[i, k] = normalized metric k for candidate i."""
    sub = df[df[f"{pname}_ok"] == True].copy()  # noqa: E712
    sub = sub[sub[f"score_{pname}"].notna()]
    if sub.empty:
        return None, None
    n = normalized_metrics(sub, pname)          # SAME code the screen scores with
    M = np.column_stack([np.asarray(n[k], float) for k in METRICS])
    return sub["name"].to_numpy(), np.nan_to_num(M, nan=0.0)


def nominal_weights(pname):
    w = PROFILES[pname]["weights"]
    v = np.array([w.get(k, 0.0) for k in METRICS], float)
    return v / v.sum()


def family_of(names):
    return pd.Series(names).str.replace(r"_v\d+$", "", regex=True)


def spearman(a, b):
    """Rank correlation, without pulling in scipy.

    NOTE: pd.Series.rank().to_numpy() can hand back a READ-ONLY view of the
    underlying block, depending on the pandas version. Mutating it in place
    (`ra -= ra.mean()`) therefore works on some machines and raises
    "ValueError: output array is read-only" on others. Build new arrays instead
    of mutating - never assume a numpy view from pandas is writable.
    """
    ra = np.asarray(pd.Series(a).rank(), dtype=float)
    rb = np.asarray(pd.Series(b).rank(), dtype=float)
    ra = ra - ra.mean()
    rb = rb - rb.mean()
    d = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / d) if d else np.nan


def analyze(names, M, w_nom, draws, conc, rng):
    s_nom = M @ w_nom
    order = np.argsort(-s_nom)
    winner = names[order[0]]
    winner_fam = family_of([winner])[0]
    margin = float(s_nom[order[0]] - s_nom[order[1]]) if len(s_nom) > 1 else np.nan

    W = rng.dirichlet(w_nom * conc, size=draws)          # (draws, n_metrics)
    S = M @ W.T                                          # (n_cand, draws)
    win_idx = S.argmax(axis=0)
    win_prob = pd.Series(names[win_idx]).value_counts() / draws
    fam_prob = family_of(names[win_idx]).value_counts() / draws

    k = min(len(names), 200)
    top = order[:k]
    rhos = np.array([spearman(s_nom[top], S[top, j])
                     for j in range(min(draws, 400))])

    oat = {}
    for m in range(len(METRICS)):
        ranks = []
        for wm in np.linspace(0.0, 1.0, 21):
            w = w_nom.copy()
            other = w_nom.sum() - w_nom[m]
            w = w * ((1.0 - wm) / other) if other > 1e-9 else np.full_like(w, (1.0 - wm) / (len(w) - 1))
            w[m] = wm
            s = M @ w
            ranks.append(int((s > s[order[0]]).sum()) + 1)
        oat[METRICS[m]] = np.array(ranks)

    vertices = {}
    for m, met in enumerate(METRICS):
        w = np.zeros(len(METRICS)); w[m] = 1.0
        vertices[met] = names[int(np.argmax(M @ w))]

    return dict(winner=winner, winner_fam=winner_fam, margin=margin, s_nom=s_nom,
                order=order, win_prob=win_prob, fam_prob=fam_prob, rhos=rhos,
                oat=oat, vertices=vertices)


def figure(pname, names, res, plotdir, draws, conc):
    w_nom = nominal_weights(pname)
    fam_p = res["fam_prob"].get(res["winner_fam"], 0.0) * 100
    var_p = res["win_prob"].iloc[0] * 100

    fig, ((a1, a2), (a3, a4)) = plt.subplots(2, 2, figsize=(12.5, 8.4))

    # (a) win probability by section
    wp = res["win_prob"].head(8)[::-1]
    cols = [C_WINNER if n == res["winner"] else C_RIVAL for n in wp.index]
    a1.barh(range(len(wp)), wp.values * 100, color=cols, height=0.72)
    a1.set_yticks(range(len(wp)))
    a1.set_yticklabels(wp.index)
    a1.set_xlabel("Win probability [%]")
    a1.set_title("(a)  Winning section across sampled weight vectors", loc="left")
    a1.grid(axis="x", color=C_GRID, lw=0.6)
    a1.set_axisbelow(True)
    a1.text(0.98, 0.06,
            f"Variant {res['winner']}:  {var_p:.0f}%\n"
            f"Family {res['winner_fam']}:  {fam_p:.0f}%",
            transform=a1.transAxes, ha="right", va="bottom", fontsize=8.5,
            bbox=dict(fc="white", ec="0.7", lw=0.7, boxstyle="round,pad=0.4"))

    # (b) one-at-a-time weight sweep
    x = np.linspace(0, 1, 21)
    for m, met in enumerate(METRICS):
        a2.plot(x, res["oat"][met], lw=1.4, marker="o", ms=2.6, label=met)
    a2.plot(w_nom, np.ones_like(w_nom), ls="none", marker="D", ms=5,
            mfc="none", mec="k", mew=1.1, label="Nominal weight")
    a2.set_yscale("log")
    a2.invert_yaxis()
    a2.set_xlabel("Weight assigned to single metric [-]")
    a2.set_ylabel(f"Rank of {res['winner']}")
    a2.set_title("(b)  Rank of nominal winner under one-at-a-time re-weighting",
                 loc="left")
    a2.legend(ncol=2, loc="lower left")
    a2.grid(color=C_GRID, lw=0.6)
    a2.set_axisbelow(True)

    # (c) nominal ranking and margin
    s = res["s_nom"][res["order"]][:15][::-1]
    nm = names[res["order"]][:15][::-1]
    cols = [C_WINNER if n == res["winner"] else C_GHOST for n in nm]
    a3.barh(range(len(s)), s, color=cols, height=0.72)
    a3.set_yticks(range(len(s)))
    a3.set_yticklabels(nm, fontsize=7.5)
    a3.set_xlim(min(s) * 0.985, max(s) * 1.004)
    a3.set_xlabel(f"score_{pname} [-]")
    a3.set_title(f"(c)  Nominal ranking, top 15   "
                 f"(margin over 2nd = {res['margin']:.4f})", loc="left")
    a3.grid(axis="x", color=C_GRID, lw=0.6)
    a3.set_axisbelow(True)

    # (d) rank-order stability
    a4.hist(res["rhos"], bins=25, color=C_RIVAL, edgecolor="white", lw=0.4)
    med = float(np.median(res["rhos"]))
    a4.axvline(med, color=C_WINNER, lw=1.6,
               label=rf"median $\rho$ = {med:.3f}")
    a4.legend(loc="upper left")
    a4.set_xlabel(r"Spearman $\rho$ vs. nominal ranking [-]")
    a4.set_ylabel("Draws")
    a4.set_title(r"(d)  Rank-order stability, top 200 candidates", loc="left")
    a4.grid(color=C_GRID, lw=0.6)
    a4.set_axisbelow(True)

    fig.suptitle(f"{pname}  —  scoring weight sensitivity", x=0.012, ha="left",
                 fontsize=13, fontweight="bold")
    fig.text(0.012, 0.945,
             f"Dirichlet(w$_{{nom}}$·κ), κ = {conc:g}, N = {draws} draws     |     "
             f"nominal winner: {res['winner']} (family {res['winner_fam']})     |     "
             f"family win probability: {fam_p:.0f}%",
             ha="left", fontsize=8.5, color=C_TEXT)
    fig.text(0.012, 0.012,
             "Single-metric optima (weight-simplex vertices):   "
             + "    ".join(f"{k} → {v}" for k, v in res["vertices"].items()),
             ha="left", fontsize=7.5, color=C_TEXT)

    fig.tight_layout(rect=[0, 0.03, 1, 0.925])
    out = os.path.join(plotdir, f"weight_sensitivity_{pname}.png")
    fig.savefig(out)
    plt.close(fig)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="screen_results")
    ap.add_argument("--draws", type=int, default=4000)
    ap.add_argument("--concentration", type=float, default=40.0,
                    help="Dirichlet concentration: high = weights stay near "
                         "nominal, low = wide disagreement (default 40)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    load_run_config(args.results)
    plotdir = os.path.join(args.results, "plots")
    os.makedirs(plotdir, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    df = pd.read_excel(os.path.join(args.results, "results.xlsx"),
                       sheet_name="all_survivors")

    print(f"Weight sensitivity: Dirichlet(w_nom * {args.concentration:g}), "
          f"{args.draws} draws\n")
    hdr = (f"{'element':<15s}{'winner':<14s}{'family':<10s}"
           f"{'variant':>8s}{'family':>8s}{'margin':>9s}{'rho':>7s}   verdict")
    print(hdr); print("-" * len(hdr))

    rows = []
    for pname in PROFILES:
        names, M = build_matrix(df, pname)
        if names is None:
            print(f"{pname:<15s} no qualifiers")
            continue
        res = analyze(names, M, nominal_weights(pname), args.draws,
                      args.concentration, rng)
        figure(pname, names, res, plotdir, args.draws, args.concentration)

        fam_p = res["fam_prob"].get(res["winner_fam"], 0.0) * 100
        var_p = res["win_prob"].iloc[0] * 100
        # Verdict is judged on the FAMILY: a mould is cut per family, not variant.
        verdict = ("robust" if fam_p >= 75 else
                   "fragile" if fam_p < 50 else "moderate")
        print(f"{pname:<15s}{res['winner']:<14s}{res['winner_fam']:<10s}"
              f"{var_p:7.0f}%{fam_p:7.0f}%{res['margin']:9.4f}"
              f"{np.median(res['rhos']):7.3f}   {verdict}")

        for name, p in res["win_prob"].head(10).items():
            rows.append(dict(profile=pname, level="variant", name=name,
                             win_prob=p, is_nominal_winner=(name == res["winner"]),
                             margin_over_2nd=res["margin"],
                             median_spearman=float(np.median(res["rhos"]))))
        for name, p in res["fam_prob"].head(6).items():
            rows.append(dict(profile=pname, level="family", name=name,
                             win_prob=p,
                             is_nominal_winner=(name == res["winner_fam"]),
                             margin_over_2nd=res["margin"],
                             median_spearman=float(np.median(res["rhos"]))))

    csv = os.path.join(args.results, "weight_sensitivity.csv")
    pd.DataFrame(rows).to_csv(csv, index=False)
    print(f"\nFigures -> {plotdir}/weight_sensitivity_<element>.png")
    print(f"Data    -> {csv}")
    print("\nVerdict is judged on FAMILY win probability: a mould is cut per")
    print("family, not per CST variant. Rivals that are siblings of the winner")
    print("do not represent a different aerodynamic decision.")


if __name__ == "__main__":
    main()
