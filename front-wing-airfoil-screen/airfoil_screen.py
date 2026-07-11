#!/usr/bin/env python3
"""
FSAE front-wing airfoil screening pipeline (NeuralFoil).

Screens the UIUC database (~2,175 airfoils, bundled with AeroSandbox) plus
CST/Kulfan-perturbed variants of the best seeds, at FSAE front-wing Reynolds
numbers, and ranks candidates against TWO spanwise role profiles:

  CENTER   - feeds the undertray: moderate loading, high L/D at the target CL
             (clean/thin wake), gentle stall, low |Cm|.
  OUTBOARD - drives multi-element downforce + outwash around the front tires:
             high CLmax, high usable CL with stall margin, gentle stall.

Outputs (in --out dir):
  results.xlsx          all metrics + ranked sheets per role + config
  plots/*.png           polars for the top-N of each role
  shortlist_dat/*.dat   Selig-format coordinates for XFLR5 validation

Usage:
  pip install neuralfoil aerosandbox pandas openpyxl matplotlib tqdm
  python airfoil_screen.py                 # full screen (~10-20 min)
  python airfoil_screen.py --quick         # 150-airfoil smoke test
  python airfoil_screen.py --no-variants   # UIUC pool only

Notes:
  - NeuralFoil is a surrogate for XFoil: treat this as a FILTER, not a final
    answer. Validate the shortlist in XFLR5/XFoil, then CFD with all elements
    + ground effect. Single-element polars can't capture slot interactions.
  - analysis_confidence < CONF_MIN points are discarded; airfoils NeuralFoil
    is unsure about are gated out rather than trusted.
"""

import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ----------------------------------------------------------------------------
# CONFIG - tune these to your car
# ----------------------------------------------------------------------------
RE_LIST = [200_000, 350_000, 500_000]   # chord Re sweep (edit for your chords/speeds)
ALPHA = np.arange(-2.0, 18.1, 0.5)      # AoA sweep, deg
MODEL_SIZE = "large"                    # NeuralFoil model: xsmall..xxxlarge

CONF_MIN = 0.85                         # min NeuralFoil analysis confidence
TC_MIN, TC_MAX = 0.055, 0.16            # thickness gates (manufacturability / drag)
CLMAX_GATE = 1.40                       # hard gate: worst-case CLmax across Re
MIN_VALID_FRAC = 0.6                    # min fraction of alpha points passing conf

CL_TARGET_CENTER = 1.40                 # center-section operating CL (undertray feed)
CL_TARGET_OUTBOARD = 1.70               # outboard operating CL

# Scoring weights (normalized 0-1 metrics; sums need not equal 1)
WEIGHTS_CENTER = {
    "LD_at_target": 0.35,   # L/D at CL_TARGET_CENTER -> thin, low-loss wake to undertray
    "CL_usable":    0.15,   # lift capability with 2 deg stall margin
    "stall_gentle": 0.25,   # soft post-stall CL drop -> robust to ride-height/yaw upwash changes
    "Cm_low":       0.15,   # low |Cm| -> less pitch sensitivity, helps aero balance
    "LD_max":       0.10,
}
WEIGHTS_OUTBOARD = {
    "CLmax":        0.35,   # peak sectional lift for the multi-element stack
    "CL_usable":    0.25,   # lift at (stall - 2 deg): what you can actually run
    "stall_gentle": 0.20,   # tolerant to flap-induced loading / tire wake dirt
    "LD_at_target": 0.10,   # at CL_TARGET_OUTBOARD
    "Cm_low":       0.10,
}

N_SEEDS_FOR_VARIANTS = 15               # top seeds (by combined score) to perturb
N_VARIANTS_PER_SEED = 30
VARIANT_SIGMA = 0.08                    # relative Gaussian noise on Kulfan weights
TOP_N = 15                              # shortlist size per role
# ----------------------------------------------------------------------------


def get_uiuc_names():
    import aerosandbox
    db = os.path.join(os.path.dirname(aerosandbox.__file__),
                      "geometry", "airfoil", "airfoil_database")
    return sorted(f[:-4] for f in os.listdir(db) if f.endswith(".dat"))


def polar_metrics(aero, alpha):
    """Extract screening metrics from one NeuralFoil polar. Returns None if unusable."""
    conf = np.asarray(aero["analysis_confidence"], float)
    CL = np.asarray(aero["CL"], float)
    CD = np.asarray(aero["CD"], float)
    CM = np.asarray(aero["CM"], float)

    valid = conf >= CONF_MIN
    if valid.mean() < MIN_VALID_FRAC:
        return None
    a, cl, cd, cm = alpha[valid], CL[valid], CD[valid], CM[valid]

    i_max = int(np.argmax(cl))
    if i_max in (0, len(cl) - 1):       # CLmax on sweep boundary -> untrustworthy
        return None
    clmax, a_stall = float(cl[i_max]), float(a[i_max])

    # usable point: 2 deg below stall
    a_use = a_stall - 2.0
    cl_use = float(np.interp(a_use, a, cl))
    cm_use = float(np.interp(a_use, a, cm))

    # stall gentleness: CL retention 3 deg past stall (0..1, 1 = no drop)
    a_post = min(a_stall + 3.0, a[-1])
    cl_post = float(np.interp(a_post, a, cl))
    gentle = max(0.0, min(1.0, cl_post / clmax))

    # pre-stall branch for CD(CL) / L/D(CL) interpolation
    pre = slice(0, i_max + 1)
    cl_pre, cd_pre = cl[pre], cd[pre]
    order = np.argsort(cl_pre)
    cl_s, cd_s = cl_pre[order], cd_pre[order]

    def ld_at(cl_t):
        if clmax < cl_t:
            return 0.0                  # can't reach target -> zero credit
        cd_t = float(np.interp(cl_t, cl_s, cd_s))
        return cl_t / max(cd_t, 1e-6)

    return dict(
        CLmax=clmax, alpha_stall=a_stall, CL_usable=cl_use,
        stall_gentle=gentle, Cm_use=cm_use,
        LD_max=float(np.max(cl_pre[1:] / np.maximum(cd_pre[1:], 1e-6))) if i_max > 1 else 0.0,
        LD_at_center=ld_at(CL_TARGET_CENTER),
        LD_at_outboard=ld_at(CL_TARGET_OUTBOARD),
    )


def evaluate(name, kulfan_params, tc, camber):
    """Run NeuralFoil at all Re, aggregate metrics. Returns row dict or None."""
    import neuralfoil as nf
    per_re = []
    for Re in RE_LIST:
        aero = nf.get_aero_from_kulfan_parameters(
            kulfan_parameters=kulfan_params, alpha=ALPHA, Re=Re,
            model_size=MODEL_SIZE)
        m = polar_metrics(aero, ALPHA)
        if m is None:
            return None                 # unusable at any Re -> gate out
        per_re.append(m)

    def mean(k): return float(np.mean([m[k] for m in per_re]))
    def worst(k): return float(np.min([m[k] for m in per_re]))

    return dict(
        name=name, t_c=tc, camber=camber,
        CLmax_worst=worst("CLmax"), CLmax_mean=mean("CLmax"),
        CL_usable=mean("CL_usable"), alpha_stall=mean("alpha_stall"),
        stall_gentle=mean("stall_gentle"), Cm_use=mean("Cm_use"),
        LD_max=mean("LD_max"),
        LD_at_CL_center=worst("LD_at_center"),
        LD_at_CL_outboard=worst("LD_at_outboard"),
    )


def kulfan_dict(kf):
    return dict(upper_weights=np.array(kf.upper_weights),
                lower_weights=np.array(kf.lower_weights),
                leading_edge_weight=float(kf.leading_edge_weight),
                TE_thickness=float(kf.TE_thickness))


def score(df):
    """Add normalized metrics + role scores. Normalization: 5th-95th pct min-max."""
    def norm(s, invert=False):
        lo, hi = s.quantile(0.05), s.quantile(0.95)
        x = ((s - lo) / max(hi - lo, 1e-9)).clip(0, 1)
        return 1 - x if invert else x

    n = dict(
        CLmax=norm(df.CLmax_worst),
        CL_usable=norm(df.CL_usable),
        stall_gentle=norm(df.stall_gentle),
        Cm_low=norm(df.Cm_use.abs(), invert=True),
        LD_max=norm(df.LD_max),
    )
    nc = dict(n); nc["LD_at_target"] = norm(df.LD_at_CL_center)
    no = dict(n); no["LD_at_target"] = norm(df.LD_at_CL_outboard)

    df["score_center"] = sum(w * nc[k] for k, w in WEIGHTS_CENTER.items()) / sum(WEIGHTS_CENTER.values())
    df["score_outboard"] = sum(w * no[k] for k, w in WEIGHTS_OUTBOARD.items()) / sum(WEIGHTS_OUTBOARD.values())
    return df


def run_pool(pool, label):
    """pool: list of (name, kulfan_params, t_c, camber). Returns DataFrame."""
    from tqdm import tqdm
    rows = []
    for name, kp, tc, camber in tqdm(pool, desc=label, ncols=80):
        try:
            r = evaluate(name, kp, tc, camber)
        except Exception:
            r = None
        if r is not None and r["CLmax_worst"] >= CLMAX_GATE:
            rows.append(r)
    return pd.DataFrame(rows)


def build_uiuc_pool(limit=None):
    import aerosandbox as asb
    names = get_uiuc_names()
    if limit:
        names = names[:limit]
    pool = []
    print(f"Preparing {len(names)} UIUC airfoils (geometry gates: "
          f"{TC_MIN} <= t/c <= {TC_MAX})...")
    for name in names:
        try:
            af = asb.Airfoil(name)
            if af.coordinates is None or len(af.coordinates) < 20:
                continue
            tc = float(af.max_thickness())
            if not (TC_MIN <= tc <= TC_MAX):
                continue
            camber = float(af.max_camber())
            if camber < 0.015:          # front wing: symmetric/low-camber pointless
                continue
            kf = af.to_kulfan_airfoil()
            pool.append((name, kulfan_dict(kf), tc, camber))
        except Exception:
            continue
    print(f"  {len(pool)} passed geometry gates.")
    return pool


def build_variant_pool(df, uiuc_pool, rng):
    import aerosandbox as asb
    lookup = {name: kp for name, kp, _, _ in uiuc_pool}
    df = df.copy()
    df["combo"] = df.score_center + df.score_outboard
    seeds = df.nlargest(N_SEEDS_FOR_VARIANTS, "combo")["name"].tolist()
    pool = []
    for seed in seeds:
        kp = lookup.get(seed)
        if kp is None:
            continue
        for j in range(N_VARIANTS_PER_SEED):
            up = kp["upper_weights"] * (1 + rng.normal(0, VARIANT_SIGMA, len(kp["upper_weights"])))
            lo = kp["lower_weights"] * (1 + rng.normal(0, VARIANT_SIGMA, len(kp["lower_weights"])))
            newkp = dict(upper_weights=up, lower_weights=lo,
                         leading_edge_weight=kp["leading_edge_weight"] * (1 + rng.normal(0, VARIANT_SIGMA)),
                         TE_thickness=kp["TE_thickness"])
            try:
                kf = asb.KulfanAirfoil(name=f"{seed}_v{j:02d}", **newkp)
                tc = float(kf.max_thickness())
                if not (TC_MIN <= tc <= TC_MAX):
                    continue
                camber = float(kf.max_camber())
                pool.append((f"{seed}_v{j:02d}", newkp, tc, camber))
            except Exception:
                continue
    print(f"Generated {len(pool)} valid variants from {len(seeds)} seeds.")
    return pool


def export_outputs(df, uiuc_pool, variant_pool, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import aerosandbox as asb
    import neuralfoil as nf

    os.makedirs(outdir, exist_ok=True)
    plots = os.path.join(outdir, "plots"); os.makedirs(plots, exist_ok=True)
    datd = os.path.join(outdir, "shortlist_dat"); os.makedirs(datd, exist_ok=True)

    top_c = df.nlargest(TOP_N, "score_center")
    top_o = df.nlargest(TOP_N, "score_outboard")

    # spreadsheet
    cols = ["name", "t_c", "camber", "CLmax_worst", "CLmax_mean", "CL_usable",
            "alpha_stall", "stall_gentle", "Cm_use", "LD_max",
            "LD_at_CL_center", "LD_at_CL_outboard", "score_center", "score_outboard"]
    cfg = pd.DataFrame(
        [("Re_list", RE_LIST), ("alpha", f"{ALPHA[0]}..{ALPHA[-1]}"),
         ("model_size", MODEL_SIZE), ("conf_min", CONF_MIN),
         ("t/c gates", (TC_MIN, TC_MAX)), ("CLmax gate", CLMAX_GATE),
         ("CL target center", CL_TARGET_CENTER),
         ("CL target outboard", CL_TARGET_OUTBOARD),
         ("weights center", WEIGHTS_CENTER), ("weights outboard", WEIGHTS_OUTBOARD)],
        columns=["parameter", "value"]).astype(str)
    xlsx = os.path.join(outdir, "results.xlsx")
    with pd.ExcelWriter(xlsx, engine="openpyxl") as w:
        df.sort_values("score_outboard", ascending=False)[cols].to_excel(w, "all_survivors", index=False)
        top_c[cols].to_excel(w, "center_top", index=False)
        top_o[cols].to_excel(w, "outboard_top", index=False)
        cfg.to_excel(w, "config", index=False)

    # polar plots + .dat export for shortlist
    lookup = {name: kp for name, kp, _, _ in (uiuc_pool + variant_pool)}
    shortlist = pd.concat([top_c, top_o]).drop_duplicates("name")
    for _, row in shortlist.iterrows():
        kp = lookup[row["name"]]
        kf = asb.KulfanAirfoil(name=row["name"], **kp)
        # Selig .dat for XFLR5
        c = kf.coordinates
        with open(os.path.join(datd, f"{row['name']}.dat"), "w") as f:
            f.write(row["name"] + "\n")
            for x, y in c:
                f.write(f" {x:.6f}  {y:.6f}\n")
        # polar plot
        fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
        for Re in RE_LIST:
            aero = nf.get_aero_from_kulfan_parameters(
                kulfan_parameters=kp, alpha=ALPHA, Re=Re, model_size=MODEL_SIZE)
            m = np.asarray(aero["analysis_confidence"]) >= CONF_MIN
            ax[0].plot(ALPHA[m], np.asarray(aero["CL"])[m], label=f"Re {Re//1000}k")
            ax[1].plot(np.asarray(aero["CD"])[m], np.asarray(aero["CL"])[m])
            ax[2].plot(ALPHA[m], np.asarray(aero["CL"])[m] / np.asarray(aero["CD"])[m])
        ax[0].set_xlabel("alpha"); ax[0].set_ylabel("CL"); ax[0].legend(fontsize=7)
        ax[1].set_xlabel("CD"); ax[1].set_ylabel("CL")
        ax[2].set_xlabel("alpha"); ax[2].set_ylabel("L/D")
        for a_ in ax: a_.grid(alpha=0.3)
        fig.suptitle(row["name"])
        fig.tight_layout()
        fig.savefig(os.path.join(plots, f"{row['name']}.png"), dpi=110)
        plt.close(fig)

    print(f"\nWrote {xlsx}")
    print(f"Wrote {len(shortlist)} polar plots -> {plots}")
    print(f"Wrote {len(shortlist)} .dat files  -> {datd}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--out", default="screen_results")
    ap.add_argument("--quick", action="store_true", help="150-airfoil smoke test")
    ap.add_argument("--no-variants", action="store_true")
    ap.add_argument("--model-size", default=None, help="override NeuralFoil model size")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    global MODEL_SIZE
    if args.model_size:
        MODEL_SIZE = args.model_size
    rng = np.random.default_rng(args.seed)

    uiuc_pool = build_uiuc_pool(limit=150 if args.quick else None)
    df = run_pool(uiuc_pool, "UIUC screen")
    if df.empty:
        sys.exit("No airfoils survived the gates - loosen CONFIG values.")
    df = score(df)
    print(f"UIUC survivors: {len(df)} / {len(uiuc_pool)}")

    variant_pool = []
    if not args.no_variants:
        variant_pool = build_variant_pool(df, uiuc_pool, rng)
        dfv = run_pool(variant_pool, "Variant screen")
        if not dfv.empty:
            df = score(pd.concat([df.drop(columns=["score_center", "score_outboard", "combo"],
                                          errors="ignore"), dfv], ignore_index=True))

    export_outputs(df, uiuc_pool, variant_pool, args.out)

    for role in ("center", "outboard"):
        print(f"\nTop 5 - {role.upper()}:")
        t = df.nlargest(5, f"score_{role}")
        for _, r in t.iterrows():
            print(f"  {r['name']:<22} CLmax={r.CLmax_worst:.2f} "
                  f"L/D@CL={r.LD_at_CL_center if role=='center' else r.LD_at_CL_outboard:6.1f} "
                  f"gentle={r.stall_gentle:.2f} Cm={r.Cm_use:+.3f} "
                  f"score={r[f'score_{role}']:.3f}")


if __name__ == "__main__":
    main()
