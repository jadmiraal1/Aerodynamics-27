#!/usr/bin/env python3
"""
FSAE front-wing airfoil screening pipeline (NeuralFoil) - per-element edition.

Screens the UIUC database (~2,175 airfoils, bundled with AeroSandbox) plus
CST/Kulfan-perturbed variants of the best seeds, and ranks candidates
separately for each wing element AND span station:

  MAIN-CENTER   (c~275mm, Re 300-510k) - center span, main plane only:
        feeds the undertray. Lightly loaded (CL~1.1), scored almost entirely
        on wake quality (L/D at target), stall robustness, and low |Cm|.
  MAIN-OUTBOARD (c~275mm, Re 300-510k) - sets the wake and the car's aero
        balance: high L/D at moderate CL, gentle stall, low |Cm|,
        thick enough for a spar.
  FLAP1 (c~170mm, Re 185-315k) - loaded by the main element's circulation:
        high usable CL, camber welcome, Cm penalty relaxed.
  FLAP2 (c~110mm, Re 120-205k) - smallest chord, lowest Re, most aggressive
        loading: CLmax-dominated, thin sections fine.

Re bands come from the TR26 lap sim (aero-weighted speeds 16-27 m/s, design
point 20 m/s) and chords from the front-bulkhead packaging study (Jul 2026).

Outputs (in --out dir):
  results.xlsx          all survivors + ranked sheet per element + config
  plots/*.png           polars for the top-N of each element
  shortlist_dat/*.dat   Selig-format coordinates for XFLR5 validation

Usage:
  pip install neuralfoil aerosandbox pandas openpyxl matplotlib tqdm
  python airfoil_screen.py                 # full screen (~1-2 min)
  python airfoil_screen.py --quick         # 150-airfoil smoke test
  python airfoil_screen.py --no-variants   # UIUC pool only

Notes:
  - NeuralFoil is a surrogate for XFoil: treat this as a FILTER, not a final
    answer. Validate the shortlist in XFLR5/XFoil, then CFD with all elements
    + ground effect. Single-element polars can't capture slot interactions -
    flap candidates especially operate in the main element's downwash field,
    which this screen approximates only through role-appropriate weighting.
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
NU = 1.46e-5                            # kinematic viscosity, sea level [m^2/s]
V_DESIGN = (16.0, 20.0, 27.0)           # lo / design / hi speeds [m/s] (lap sim)

# Element profiles. re_list is computed from chord x speeds; edit chords here.
PROFILES = {
    "main_center": dict(
        chord=0.275,
        cl_target=1.10,                 # lightly loaded: protect the undertray inlet
        clmax_gate=1.30,                # still want stall margin above CL 1.1
        tc=(0.080, 0.160),              # same spar carries through
        weights={
            "LD_at_target": 0.40,       # wake thinness at the undertray-feed loading
            "stall_gentle": 0.25,       # closest element to ground at the inlet
            "Cm_low":       0.20,       # balance + ride-height pitch sensitivity
            "CL_usable":    0.05,
            "LD_max":       0.10,
        },
    ),
    "main_outboard": dict(
        chord=0.275,
        cl_target=1.50,                 # operating CL for L/D scoring
        clmax_gate=1.40,                # hard gate on worst-case CLmax
        tc=(0.080, 0.160),              # spar needs thickness
        weights={
            "LD_at_target": 0.35,       # thin, low-loss wake -> undertray + balance
            "CL_usable":    0.15,
            "stall_gentle": 0.25,       # robust to ride height / pitch / yaw
            "Cm_low":       0.15,       # pitch sensitivity, structural twist
            "LD_max":       0.10,
        },
    ),
    "flap1": dict(
        chord=0.170,
        cl_target=1.80,
        clmax_gate=1.60,
        tc=(0.055, 0.130),
        weights={
            "CLmax":        0.35,
            "CL_usable":    0.25,       # lift at 2 deg stall margin
            "stall_gentle": 0.20,       # tolerant to main-element upwash changes
            "LD_at_target": 0.10,
            "Cm_low":       0.10,       # relaxed: flap Cm reacts through the stack
        },
    ),
    "flap2": dict(
        chord=0.110,
        cl_target=1.75,                 # was 1.90: sat ~at the qualifier CLmax band,
                                        # so L/D was interpolated at the stall knee
                                        # (only 15% of qualifiers could reach it).
                                        # 1.75 -> ~50% reach, +0.10 median margin,
                                        # matching flap1's healthy spread.
        clmax_gate=1.60,
        tc=(0.045, 0.120),
        weights={
            "CLmax":        0.40,
            "CL_usable":    0.25,
            "stall_gentle": 0.20,
            "LD_at_target": 0.10,
            "Cm_low":       0.05,
        },
    ),
}
for _p in PROFILES.values():            # Re list from chord x (lo, design, hi)
    _p["re_list"] = [int(v * _p["chord"] / NU) for v in V_DESIGN]

ALPHA = np.arange(-2.0, 18.1, 0.5)      # AoA sweep, deg
MODEL_SIZE = "large"                    # NeuralFoil model: xsmall..xxxlarge
CONF_MIN = 0.85                         # min NeuralFoil analysis confidence
MIN_VALID_FRAC = 0.6                    # min fraction of alpha points passing conf
CAMBER_MIN = 0.015                      # symmetric sections pointless on a wing

N_SEEDS_PER_PROFILE = 6                 # top seeds per element to perturb
N_VARIANTS_PER_SEED = 30
VARIANT_SIGMA = 0.08                    # relative Gaussian noise on Kulfan weights
TOP_N = 15                              # shortlist size per element
# ----------------------------------------------------------------------------

TC_MIN_ALL = min(p["tc"][0] for p in PROFILES.values())
TC_MAX_ALL = max(p["tc"][1] for p in PROFILES.values())
ALL_RE = sorted({re for p in PROFILES.values() for re in p["re_list"]})
CL_TARGETS = {name: p["cl_target"] for name, p in PROFILES.items()}


def get_uiuc_names():
    import aerosandbox
    db = os.path.join(os.path.dirname(aerosandbox.__file__),
                      "geometry", "airfoil", "airfoil_database")
    return sorted(f[:-4] for f in os.listdir(db) if f.endswith(".dat"))


def polar_metrics(aero, alpha):
    """Metrics from one NeuralFoil polar (one Re). None if unusable."""
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

    a_use = a_stall - 2.0               # usable point: 2 deg below stall
    cl_use = float(np.interp(a_use, a, cl))
    cm_use = float(np.interp(a_use, a, cm))

    a_post = min(a_stall + 3.0, a[-1])  # CL retention 3 deg past stall (0..1)
    gentle = max(0.0, min(1.0, float(np.interp(a_post, a, cl)) / clmax))

    pre = slice(0, i_max + 1)           # pre-stall branch for L/D(CL)
    cl_pre, cd_pre = cl[pre], cd[pre]
    order = np.argsort(cl_pre)
    cl_s, cd_s = cl_pre[order], cd_pre[order]

    def ld_at(cl_t):
        if clmax < cl_t:
            return 0.0                  # can't reach target -> zero credit
        return cl_t / max(float(np.interp(cl_t, cl_s, cd_s)), 1e-6)

    return dict(
        CLmax=clmax, alpha_stall=a_stall, CL_usable=cl_use,
        stall_gentle=gentle, Cm_use=cm_use,
        LD_max=float(np.max(cl_pre[1:] / np.maximum(cd_pre[1:], 1e-6))) if i_max > 1 else 0.0,
        LD_at={name: ld_at(t) for name, t in CL_TARGETS.items()},
    )


def evaluate(name, kulfan_params, tc, camber):
    """NeuralFoil at every unique Re; aggregate per element profile."""
    import neuralfoil as nf
    by_re = {}
    for Re in ALL_RE:
        aero = nf.get_aero_from_kulfan_parameters(
            kulfan_parameters=kulfan_params, alpha=ALPHA, Re=Re,
            model_size=MODEL_SIZE)
        by_re[Re] = polar_metrics(aero, ALPHA)

    row = dict(name=name, t_c=tc, camber=camber)
    any_valid = False
    for pname, p in PROFILES.items():
        ms = [by_re[Re] for Re in p["re_list"]]
        ok = (all(m is not None for m in ms)
              and p["tc"][0] <= tc <= p["tc"][1]
              and min(m["CLmax"] for m in ms) >= p["clmax_gate"])
        row[f"{pname}_ok"] = ok
        if not ok:
            continue
        any_valid = True
        row[f"{pname}_CLmax"] = min(m["CLmax"] for m in ms)          # worst case
        row[f"{pname}_CL_usable"] = float(np.mean([m["CL_usable"] for m in ms]))
        row[f"{pname}_stall_gentle"] = float(np.mean([m["stall_gentle"] for m in ms]))
        row[f"{pname}_Cm"] = float(np.mean([m["Cm_use"] for m in ms]))
        row[f"{pname}_LD_max"] = float(np.mean([m["LD_max"] for m in ms]))
        row[f"{pname}_LD_at_CL"] = min(m["LD_at"][pname] for m in ms)  # worst case
        row[f"{pname}_alpha_stall"] = float(np.mean([m["alpha_stall"] for m in ms]))
    return row if any_valid else None


def kulfan_dict(kf):
    return dict(upper_weights=np.array(kf.upper_weights),
                lower_weights=np.array(kf.lower_weights),
                leading_edge_weight=float(kf.leading_edge_weight),
                TE_thickness=float(kf.TE_thickness))


def score(df):
    """Per-profile weighted score over 5th-95th pct min-max normalized metrics."""
    for pname, p in PROFILES.items():
        sub = df[df[f"{pname}_ok"] == True]  # noqa: E712
        if sub.empty:
            df[f"score_{pname}"] = np.nan
            continue

        def norm(s, invert=False):
            lo, hi = s.quantile(0.05), s.quantile(0.95)
            x = ((s - lo) / max(hi - lo, 1e-9)).clip(0, 1)
            return 1 - x if invert else x

        n = {
            "CLmax":        norm(sub[f"{pname}_CLmax"]),
            "CL_usable":    norm(sub[f"{pname}_CL_usable"]),
            "stall_gentle": norm(sub[f"{pname}_stall_gentle"]),
            "Cm_low":       norm(sub[f"{pname}_Cm"].abs(), invert=True),
            "LD_max":       norm(sub[f"{pname}_LD_max"]),
            "LD_at_target": norm(sub[f"{pname}_LD_at_CL"]),
        }
        w = p["weights"]
        df.loc[sub.index, f"score_{pname}"] = (
            sum(wt * n[k] for k, wt in w.items()) / sum(w.values()))
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
        if r is not None:
            rows.append(r)
    return pd.DataFrame(rows)


def build_uiuc_pool(limit=None):
    import aerosandbox as asb
    names = get_uiuc_names()
    if limit:
        names = names[:limit]
    pool = []
    print(f"Preparing {len(names)} UIUC airfoils (geometry gates: "
          f"{TC_MIN_ALL} <= t/c <= {TC_MAX_ALL}, camber >= {CAMBER_MIN})...")
    for name in names:
        try:
            af = asb.Airfoil(name)
            if af.coordinates is None or len(af.coordinates) < 20:
                continue
            tc = float(af.max_thickness())
            if not (TC_MIN_ALL <= tc <= TC_MAX_ALL):
                continue
            camber = float(af.max_camber())
            if camber < CAMBER_MIN:
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
    seeds = []
    for pname in PROFILES:
        col = f"score_{pname}"
        if col in df:
            seeds += df.nlargest(N_SEEDS_PER_PROFILE, col)["name"].tolist()
    seeds = list(dict.fromkeys(seeds))   # dedupe, keep order
    pool = []
    for seed in seeds:
        kp = lookup.get(seed)
        if kp is None:
            continue
        for j in range(N_VARIANTS_PER_SEED):
            newkp = dict(
                upper_weights=kp["upper_weights"] * (1 + rng.normal(0, VARIANT_SIGMA, len(kp["upper_weights"]))),
                lower_weights=kp["lower_weights"] * (1 + rng.normal(0, VARIANT_SIGMA, len(kp["lower_weights"]))),
                leading_edge_weight=kp["leading_edge_weight"] * (1 + rng.normal(0, VARIANT_SIGMA)),
                TE_thickness=kp["TE_thickness"])
            try:
                kf = asb.KulfanAirfoil(name=f"{seed}_v{j:02d}", **newkp)
                tc = float(kf.max_thickness())
                if not (TC_MIN_ALL <= tc <= TC_MAX_ALL):
                    continue
                pool.append((f"{seed}_v{j:02d}", newkp, tc, float(kf.max_camber())))
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

    tops = {p: df.nlargest(TOP_N, f"score_{p}") for p in PROFILES if f"score_{p}" in df}

    base = ["name", "t_c", "camber"]
    cfg_rows = [("V lo/design/hi [m/s]", V_DESIGN), ("alpha", f"{ALPHA[0]}..{ALPHA[-1]}"),
                ("model_size", MODEL_SIZE), ("conf_min", CONF_MIN)]
    for pname, p in PROFILES.items():
        cfg_rows += [(f"{pname}: chord/Re", (p["chord"], p["re_list"])),
                     (f"{pname}: CL target / CLmax gate / t/c", (p["cl_target"], p["clmax_gate"], p["tc"])),
                     (f"{pname}: weights", p["weights"])]
    cfg = pd.DataFrame(cfg_rows, columns=["parameter", "value"]).astype(str)

    xlsx = os.path.join(outdir, "results.xlsx")
    with pd.ExcelWriter(xlsx, engine="openpyxl") as w:
        df.to_excel(w, sheet_name="all_survivors", index=False)
        for pname, top in tops.items():
            cols = base + [c for c in df.columns if c.startswith(pname)] + [f"score_{pname}"]
            top[cols].to_excel(w, sheet_name=f"{pname}_top", index=False)
        cfg.to_excel(w, sheet_name="config", index=False)

    lookup = {name: kp for name, kp, _, _ in (uiuc_pool + variant_pool)}
    done = set()
    for pname, top in tops.items():
        re_list = PROFILES[pname]["re_list"]
        for _, row in top.iterrows():
            tag = f"{row['name']}__{pname}"
            if row["name"] not in done:
                kf = asb.KulfanAirfoil(name=row["name"], **lookup[row["name"]])
                with open(os.path.join(datd, f"{row['name']}.dat"), "w") as f:
                    f.write(row["name"] + "\n")
                    for x, y in kf.coordinates:
                        f.write(f" {x:.6f}  {y:.6f}\n")
                done.add(row["name"])
            fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
            for Re in re_list:
                aero = nf.get_aero_from_kulfan_parameters(
                    kulfan_parameters=lookup[row["name"]], alpha=ALPHA, Re=Re,
                    model_size=MODEL_SIZE)
                m = np.asarray(aero["analysis_confidence"]) >= CONF_MIN
                ax[0].plot(ALPHA[m], np.asarray(aero["CL"])[m], label=f"Re {Re//1000}k")
                ax[1].plot(np.asarray(aero["CD"])[m], np.asarray(aero["CL"])[m])
                ax[2].plot(ALPHA[m], np.asarray(aero["CL"])[m] / np.asarray(aero["CD"])[m])
            ax[0].set_xlabel("alpha"); ax[0].set_ylabel("CL"); ax[0].legend(fontsize=7)
            ax[1].set_xlabel("CD"); ax[1].set_ylabel("CL")
            ax[2].set_xlabel("alpha"); ax[2].set_ylabel("L/D")
            for a_ in ax: a_.grid(alpha=0.3)
            fig.suptitle(f"{row['name']}  [{pname}]")
            fig.tight_layout()
            fig.savefig(os.path.join(plots, f"{tag}.png"), dpi=110)
            plt.close(fig)

    print(f"\nWrote {xlsx}")
    print(f"Wrote polar plots -> {plots}")
    print(f"Wrote {len(done)} .dat files -> {datd}")


def main():
    ap = argparse.ArgumentParser(description="Per-element FSAE front-wing airfoil screen")
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

    print("Element profiles:")
    for pname, p in PROFILES.items():
        print(f"  {pname:6s} c={p['chord']*1000:.0f}mm  Re={p['re_list']}  "
              f"CL_target={p['cl_target']}  gate={p['clmax_gate']}  t/c={p['tc']}")

    uiuc_pool = build_uiuc_pool(limit=150 if args.quick else None)
    df = run_pool(uiuc_pool, "UIUC screen")
    if df.empty:
        sys.exit("No airfoils survived the gates - loosen CONFIG values.")
    df = score(df)
    print(f"UIUC survivors: {len(df)} / {len(uiuc_pool)}  "
          + " ".join(f"{p}:{int(df[f'{p}_ok'].sum())}" for p in PROFILES))

    variant_pool = []
    if not args.no_variants:
        variant_pool = build_variant_pool(df, uiuc_pool, rng)
        dfv = run_pool(variant_pool, "Variant screen")
        if not dfv.empty:
            df = score(pd.concat(
                [df.drop(columns=[c for c in df.columns if c.startswith("score_")]),
                 dfv], ignore_index=True))

    export_outputs(df, uiuc_pool, variant_pool, args.out)

    for pname in PROFILES:
        col = f"score_{pname}"
        if col not in df:
            continue
        print(f"\nTop 5 - {pname.upper()} (c={PROFILES[pname]['chord']*1000:.0f}mm):")
        for _, r in df.nlargest(5, col).iterrows():
            print(f"  {r['name']:<22} CLmax={r[f'{pname}_CLmax']:.2f} "
                  f"L/D@CL{PROFILES[pname]['cl_target']}={r[f'{pname}_LD_at_CL']:6.1f} "
                  f"gentle={r[f'{pname}_stall_gentle']:.2f} Cm={r[f'{pname}_Cm']:+.3f} "
                  f"score={r[col]:.3f}")


if __name__ == "__main__":
    main()
