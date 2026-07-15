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

THE ABUSE GATE
--------------
Every candidate is additionally analyzed with the leading edge tripped
(xtr=0.05 - rubber pickup, bugs, rain), at the LOW-Re end of its band: the
worst corner the car actually operates in. A section is DISQUALIFIED unless

    CLmax_abuse >= cl_target + abuse_margin

i.e. it must still reach its design loading when it gets no laminar flow.
This gate exists because the clean screen was ranking sections that its own
robustness tool knew were dead: e423 topped flap2 and main_outboard, yet
tripped it could not reach its target CL at all. A section that only works in
a clean tunnel is not a candidate, so it is removed here rather than
discovered later in a CSV nobody re-reads.

Config, element profiles, and the polar-reduction math live in common.py -
they are shared with compare.py and robustness.py so the three scripts cannot
drift apart.

Outputs (in --out dir):
  results.xlsx          all survivors + ranked sheet per element + config
  run_config.json       the config THIS run used (downstream scripts read it)
  kulfan_params.json    exact ranked geometry (downstream scripts analyze it)
  plots/*.png           polars for the top-N of each element
  shortlist_dat/*.dat   Selig-format coordinates for XFLR5 validation

Usage:
  pip install -r requirements.txt
  python airfoil_screen.py                 # full screen
  python airfoil_screen.py --quick         # 150-airfoil smoke test
  python airfoil_screen.py --no-variants   # UIUC pool only
  python airfoil_screen.py --no-abuse-gate # report abuse metrics but don't gate

Notes:
  - NeuralFoil is a surrogate for XFoil: treat this as a FILTER, not a final
    answer. Validate the shortlist in XFLR5/XFoil, then CFD with all elements
    + ground effect. Single-element polars can't capture slot interactions -
    flap candidates especially operate in the main element's downwash field,
    which this screen approximates only through role-appropriate weighting.
  - Ground effect is NOT modeled anywhere in this pipeline. For a front wing at
    h/c ~ 0.1-0.3 that is a real omission; the ranking can reorder under GE.
"""

import argparse
import json
import os
import shutil
import sys
import warnings

import numpy as np
import pandas as pd

from common import (
    AFT_X, ALL_RE, ALPHA, CL_TARGETS, CONFIG, PROFILES, TC_MAX_ALL, TC_MIN_ALL,
    V_DESIGN, abuse_kwargs, gate_kwargs, ld_at, masked_polar,
    normalized_metrics, polar_metrics, save_kulfan, save_run_config,
    truncation_station,
)


def aft_profile(kf):
    """Local thickness / chord at AFT_X. Pure geometry - computed once, at pool
    build, then reused for every element (each with its own chord)."""
    return np.array([float(kf.local_thickness(x_over_c=x)) for x in AFT_X])

warnings.filterwarnings("ignore")

N_SEEDS_PER_PROFILE = 6                 # top seeds per element to perturb
N_VARIANTS_PER_SEED = 30
VARIANT_SIGMA = 0.08                    # relative Gaussian noise on Kulfan weights
TOP_N = 15                              # shortlist size per element

# The abuse case is evaluated at the LOW-Re end of each element's band - the
# worst corner (laminar separation bubbles are nastiest at low Re, and a
# tripped LE removes the laminar run that was papering over them).
ABUSE_RE = sorted({p["re_list"][0] for p in PROFILES.values()})

APPLY_ABUSE_GATE = True                 # --no-gate flips this
APPLY_MFG_GATE = True                   # --no-mfg-gate flips this


def get_uiuc_names():
    import aerosandbox
    db = os.path.join(os.path.dirname(aerosandbox.__file__),
                      "geometry", "airfoil", "airfoil_database")
    return sorted(f[:-4] for f in os.listdir(db) if f.endswith(".dat"))


def evaluate(name, kulfan_params, tc, camber, aft=None):
    """NeuralFoil polars per element, on the BLUNTED (manufacturable) geometry.

    WHY THE GEOMETRY IS PER-ELEMENT
    -------------------------------
    These wings are hot-wire cut from foam, and a knife-edge trailing edge
    cannot be cut - the wire melts through it and the foam crumbles. Real
    airfoils have a near-zero TE, so EVERY candidate must be blunted to a
    manufacturable edge (CONFIG['te_min_mm'], default 1.5 mm) before it can be
    built.

    Crucially, "manufacturable" is a statement about ABSOLUTE thickness, so it
    depends on the CHORD - and each element has a different one. A 1.5 mm edge
    is 0.55% of main_center's 275 mm chord but 1.36% of flap2's 110 mm chord.
    The same airfoil is therefore a DIFFERENT (more blunted) section at flap2
    than at main_center, and must be analyzed separately for each. That is why
    polars are no longer shared across elements.

    WHY BLUNT RATHER THAN GATE OR TRUNCATE
    --------------------------------------
    Gating on manufacturability would rank geometry we cannot build: measured
    on the un-blunted section, then built blunted. Truncating (chopping the aft
    off) is worse still for the high-lift sections - s1223's thin cusped aft
    region IS its aft loading, and lopping 8% of chord off it destroys exactly
    the thing that made it win.

    Blunting via the Kulfan TE_thickness parameter keeps the chord and the
    camber line, adds a small base, and - because NeuralFoil evaluates
    TE_thickness natively - lets the aerodynamic COST of being manufacturable
    show up honestly in the score. Sections that tolerate a blunt edge rise;
    sections that depend on a feather edge fall. Nobody is disqualified for
    manufacturability, because every candidate is MADE manufacturable first.
    """
    import neuralfoil as nf

    def nf_call(kp, Re, **kw):
        return nf.get_aero_from_kulfan_parameters(
            kulfan_parameters=kp, alpha=ALPHA, Re=Re,
            model_size=CONFIG["model_size"], **kw)

    def reduce_case(m, target):
        """(CLmax, L/D at target) from a masked polar. NaNs if unusable."""
        if m is None:
            return np.nan, np.nan
        a_, cl_, cd_, cm_, i_ = m
        return float(cl_[i_]), ld_at(cl_, cd_, i_, target)  # NaN if unreachable

    def band_min(metrics_list, key, pname=None):
        """Worst case across the Re band; NaN if any Re is unusable."""
        if any(m is None for m in metrics_list):
            return np.nan
        vals = [(m["LD_at"][pname] if pname else m[key]) for m in metrics_list]
        return float(np.min(vals))          # NaN propagates = honest worst case

    row = dict(name=name, t_c=tc, camber=camber)
    any_valid = False
    for pname, p in PROFILES.items():
        chord_mm = p["chord"] * 1000.0
        target = p["cl_target"]
        re_lo = p["re_list"][0]

        # --- BLUNT THE SECTION TO A CUTTABLE EDGE, AT THIS ELEMENT'S CHORD ---
        te_needed = CONFIG["te_min_mm"] / chord_mm       # as a fraction of chord
        te_as_drawn = float(kulfan_params["TE_thickness"])
        te_applied = max(te_as_drawn, te_needed) if APPLY_MFG_GATE else te_as_drawn
        kp = dict(kulfan_params, TE_thickness=te_applied)

        # Polars on the BLUNTED geometry - the section we would actually cut.
        ms = [polar_metrics(nf_call(kp, Re), ALPHA) for Re in p["re_list"]]
        gms = [polar_metrics(nf_call(kp, Re, **gate_kwargs()), ALPHA)
               for Re in p["re_list"]]
        am = masked_polar(nf_call(kp, re_lo, **abuse_kwargs()), ALPHA)

        gm_lo = gms[0]
        clmax_gate_ = np.nan if gm_lo is None else gm_lo["CLmax"]
        clmax_abuse, ld_abuse = reduce_case(am, target)

        # THE GATE is the racing condition, at the low-Re end (worst corner).
        gate_ok = (not np.isnan(clmax_gate_)
                   and clmax_gate_ >= target + CONFIG["gate_margin"])

        ok = (all(m is not None for m in ms)
              and all(m is not None for m in gms)   # must be analyzable where it races
              and p["tc"][0] <= tc <= p["tc"][1]
              and min(m["CLmax"] for m in ms) >= p["clmax_gate"]
              and (gate_ok or not APPLY_ABUSE_GATE))

        row[f"{pname}_ok"] = ok
        row[f"{pname}_gate_ok"] = bool(gate_ok)

        # --- manufacturability, REPORTED (nobody is disqualified: the section
        # was MADE manufacturable above, and paid for it in the polars) ---
        row[f"{pname}_te_as_drawn_mm"] = te_as_drawn * chord_mm
        row[f"{pname}_te_applied_mm"] = te_applied * chord_mm
        # How feathered was the ORIGINAL aft? i.e. how much chord we WOULD have
        # had to chop, had we truncated instead of blunted. Pure diagnostic - a
        # high value means this section leans hard on a thin aft region, and its
        # blunted polars will differ most from the published ones.
        if aft is not None:
            x_trunc, chord_loss, _ = truncation_station(aft, p["chord"])
            row[f"{pname}_x_trunc"] = x_trunc
            row[f"{pname}_chord_loss_pct"] = 100.0 * chord_loss
        if not ok:
            continue
        any_valid = True

        clmax_clean_lo = ms[0]["CLmax"]              # same Re as the abuse case
        ld_clean_lo = ms[0]["LD_at"][pname]

        row[f"{pname}_CLmax"] = min(m["CLmax"] for m in ms)          # worst case
        row[f"{pname}_CL_usable"] = float(np.mean([m["CL_usable"] for m in ms]))
        row[f"{pname}_stall_gentle"] = float(np.mean([m["stall_gentle"] for m in ms]))
        row[f"{pname}_Cm"] = float(np.mean([m["Cm_use"] for m in ms]))
        row[f"{pname}_alpha_stall"] = float(np.mean([m["alpha_stall"] for m in ms]))

        # --- SCORED: mean L/D across the OPERATING BAND, at the race condition.
        # Not a single point: the operating CL is set by incidence (a shim at the
        # track) and inherits wide uncertainty from the lap-sim chain. Scoring at
        # one CL let that undocumented constant pick the airfoil. Worst case is
        # taken across the Re band, as everywhere else. ---
        row[f"{pname}_LD_band"] = float(np.min([m["LD_band"][pname] for m in gms]))
        row[f"{pname}_LD_flat"] = float(np.min([m["LD_flat"][pname] for m in gms]))
        row[f"{pname}_LD_max"] = float(np.mean([m["LD_max"] for m in gms]))

        # --- REPORTED: single-point L/D at the nominal target, race condition. ---
        row[f"{pname}_LD_at_CL"] = band_min(gms, None, pname)
        row[f"{pname}_LD_unreachable"] = bool(np.isnan(row[f"{pname}_LD_at_CL"]))

        # --- REPORTED: clean-case L/D, for comparison with XFLR5/published data.
        # The gap between this and *_LD_at_CL is exactly how much of a section's
        # drag advantage was laminar flow it will not have on track. ---
        row[f"{pname}_LD_at_CL_clean"] = band_min(ms, None, pname)
        row[f"{pname}_CLmax_gate"] = clmax_gate_

        # --- abuse case: reported only. A low retention here means the section
        # is fragile in the rain; it does NOT mean it was disqualified. ---
        row[f"{pname}_CLmax_abuse"] = clmax_abuse
        row[f"{pname}_LD_abuse"] = ld_abuse
        row[f"{pname}_CLmax_retention"] = (clmax_abuse / clmax_clean_lo
                                           if clmax_clean_lo else np.nan)
        row[f"{pname}_LD_retention"] = (ld_abuse / ld_clean_lo
                                        if ld_clean_lo else np.nan)
    return row if any_valid else None


def top_n(df, col, n):
    """Top n by `col`, EXCLUDING non-qualifiers.

    pandas .nlargest() pads its result with NaN rows when fewer than n non-null
    values exist. Called naively that silently injects disqualified airfoils
    into the shortlist sheets and - worse - seeds them into the variant pool.
    It never bit while every element had 200+ qualifiers; the abuse gate thins
    the flap fields enough to reach it.
    """
    return df[df[col].notna()].nlargest(n, col)


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

        # ONE definition of a normalized metric, shared with weight_sensitivity.py.
        # (Band-mean L/D already folds "unreachable" in as zero credit at each CL
        # sample inside ld_band_mean(), so no fillna() is needed here.)
        n = normalized_metrics(sub, pname)
        w = p["weights"]
        df.loc[sub.index, f"score_{pname}"] = (
            sum(wt * n[k] for k, wt in w.items()) / sum(w.values()))
    return df


def run_pool(pool, label):
    """pool: list of (name, kulfan_params, t_c, camber). Returns DataFrame."""
    from tqdm import tqdm
    rows = []
    for name, kp, tc, camber, aft in tqdm(pool, desc=label, ncols=80):
        try:
            r = evaluate(name, kp, tc, camber, aft)
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
          f"{TC_MIN_ALL} <= t/c <= {TC_MAX_ALL}, camber >= {CONFIG['camber_min']})...")
    for name in names:
        try:
            af = asb.Airfoil(name)
            if af.coordinates is None or len(af.coordinates) < 20:
                continue
            tc = float(af.max_thickness())
            if not (TC_MIN_ALL <= tc <= TC_MAX_ALL):
                continue
            camber = float(af.max_camber())
            if camber < CONFIG["camber_min"]:
                continue
            kf = af.to_kulfan_airfoil()
            pool.append((name, kulfan_dict(kf), tc, camber, aft_profile(kf)))
        except Exception:
            continue
    print(f"  {len(pool)} passed geometry gates.")
    return pool


def build_variant_pool(df, uiuc_pool, rng):
    import aerosandbox as asb
    lookup = {e[0]: e[1] for e in uiuc_pool}
    seeds = []
    for pname in PROFILES:
        col = f"score_{pname}"
        if col in df:
            seeds += top_n(df, col, N_SEEDS_PER_PROFILE)["name"].tolist()
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
                pool.append((f"{seed}_v{j:02d}", newkp, tc,
                             float(kf.max_camber()), aft_profile(kf)))
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

    # PURGE the generated dirs before writing. Without this, every run leaves
    # its artifacts behind and they accumulate: a real run of this pipeline had
    # 150 polar plots where 60 belonged, and - far worse - 102 STALE .dat files
    # from before the TE blunt existed, sitting indistinguishably beside the 60
    # current ones. Those stale files are knife-edge sections that were never
    # screened and cannot be cut. Anyone opening one in XFLR5 would validate the
    # wrong geometry and never know.
    #
    # Generated output is derived data. It should be rebuilt, not layered.
    plots = os.path.join(outdir, "plots")
    datd = os.path.join(outdir, "shortlist_dat")
    for d in (plots, datd):
        if os.path.isdir(d):
            shutil.rmtree(d)
        os.makedirs(d)
    # Figures written into plots/ by compare.py / robustness.py /
    # weight_sensitivity.py are regenerated downstream by run_pipeline.py, so
    # wiping the directory here is safe.

    lookup = {e[0]: e[1] for e in (uiuc_pool + variant_pool)}

    # Persist the run config and the EXACT ranked geometry. compare.py and
    # robustness.py consume these instead of re-deriving anything.
    save_run_config(outdir)
    save_kulfan(outdir, {n: lookup[n] for n in df["name"] if n in lookup})

    tops = {p: top_n(df, f"score_{p}", TOP_N)
            for p in PROFILES if f"score_{p}" in df}
    for p, t in tops.items():
        if len(t) < TOP_N:
            print(f"  ! {p}: only {len(t)} qualifiers (asked for {TOP_N}) - "
                  f"the gates may be too tight for this element")

    base = ["name", "t_c", "camber"]
    cfg_rows = [("V lo/design/hi [m/s]", V_DESIGN), ("alpha", f"{ALPHA[0]}..{ALPHA[-1]}"),
                ("gate applied", APPLY_ABUSE_GATE),
                ("gate case", f"n_crit={CONFIG['gate_n_crit']}, xtr={CONFIG['gate_xtr']} "
                              f"(racing condition, low-Re end)"),
                ("abuse case", f"n_crit={CONFIG['abuse_n_crit']}, xtr={CONFIG['abuse_xtr']} "
                               f"(REPORTED ONLY - never disqualifies)")]
    cfg_rows += [(k, v) for k, v in CONFIG.items()]
    for pname, p in PROFILES.items():
        cfg_rows += [(f"{pname}: chord/Re", (p["chord"], p["re_list"])),
                     (f"{pname}: CL target / CLmax gate / t/c", (p["cl_target"], p["clmax_gate"], p["tc"])),
                     (f"{pname}: gate (CLmax_gate >=)", p["cl_target"] + CONFIG["gate_margin"]),
                     (f"{pname}: weights", p["weights"])]
    cfg = pd.DataFrame(cfg_rows, columns=["parameter", "value"]).astype(str)

    xlsx = os.path.join(outdir, "results.xlsx")
    with pd.ExcelWriter(xlsx, engine="openpyxl") as w:
        df.to_excel(w, sheet_name="all_survivors", index=False)
        for pname, top in tops.items():
            cols = base + [c for c in df.columns if c.startswith(pname)] + [f"score_{pname}"]
            top[cols].to_excel(w, sheet_name=f"{pname}_top", index=False)
        cfg.to_excel(w, sheet_name="config", index=False)

    # EXPORT THE BLUNTED SECTION - the one that was screened, and the one that
    # will be cut. kulfan_params.json necessarily holds the AS-DRAWN geometry
    # (the blunt is per-element, so a single per-name file cannot carry it), and
    # exporting that raw would hand XFLR5 and the eye a knife-edge section that
    # was never analyzed: s1223_v19 at flap2's 110 mm chord is 0.13 mm at the TE
    # as drawn, and 1.50 mm as screened. Every .dat and every polar plot below
    # is therefore built from the BLUNTED geometry, and named per element,
    # because the same airfoil is a different section at 275 mm and at 110 mm.
    n_dat = 0
    for pname, top in tops.items():
        re_list = PROFILES[pname]["re_list"]
        chord_mm = PROFILES[pname]["chord"] * 1000.0
        te_needed = CONFIG["te_min_mm"] / chord_mm
        for _, row in top.iterrows():
            tag = f"{row['name']}__{pname}"
            kp = dict(lookup[row["name"]])
            te_as_drawn = float(kp["TE_thickness"])
            if APPLY_MFG_GATE:
                kp["TE_thickness"] = max(te_as_drawn, te_needed)
            te_mm = kp["TE_thickness"] * chord_mm

            kf = asb.KulfanAirfoil(name=row["name"], **kp)
            with open(os.path.join(datd, f"{tag}.dat"), "w") as f:
                f.write(f"{row['name']} [{pname}, c={chord_mm:.0f}mm, "
                        f"TE={te_mm:.2f}mm]\n")
                for x, y in kf.coordinates:
                    f.write(f" {x:.6f}  {y:.6f}\n")
            n_dat += 1

            fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
            for Re in re_list:
                aero = nf.get_aero_from_kulfan_parameters(
                    kulfan_parameters=kp, alpha=ALPHA, Re=Re,
                    model_size=CONFIG["model_size"])
                m = np.asarray(aero["analysis_confidence"]) >= CONFIG["conf_min"]
                ax[0].plot(ALPHA[m], np.asarray(aero["CL"])[m], label=f"Re {Re//1000}k")
                ax[1].plot(np.asarray(aero["CD"])[m], np.asarray(aero["CL"])[m])
                ax[2].plot(ALPHA[m], np.asarray(aero["CL"])[m] / np.asarray(aero["CD"])[m])
            ax[0].set_xlabel("alpha"); ax[0].set_ylabel("CL"); ax[0].legend(fontsize=7)
            ax[1].set_xlabel("CD"); ax[1].set_ylabel("CL")
            ax[2].set_xlabel("alpha"); ax[2].set_ylabel("L/D")
            for a_ in ax: a_.grid(alpha=0.3)
            fig.suptitle(f"{row['name']}  [{pname}]   "
                         f"c={chord_mm:.0f}mm, TE={te_mm:.2f}mm (as screened)")
            fig.tight_layout()
            fig.savefig(os.path.join(plots, f"{tag}.png"), dpi=110)
            plt.close(fig)

    print(f"\nWrote {xlsx}")
    print(f"Wrote run_config.json + kulfan_params.json -> {outdir}")
    print(f"Wrote polar plots -> {plots}")
    print(f"Wrote {n_dat} .dat files (BLUNTED, per element) -> {datd}")


def main():
    global APPLY_ABUSE_GATE, APPLY_MFG_GATE
    ap = argparse.ArgumentParser(description="Per-element FSAE front-wing airfoil screen")
    ap.add_argument("--out", default="screen_results")
    ap.add_argument("--quick", action="store_true", help="150-airfoil smoke test")
    ap.add_argument("--no-variants", action="store_true")
    ap.add_argument("--model-size", default=None, help="override NeuralFoil model size")
    ap.add_argument("--no-gate", action="store_true",
                    help="compute the gate/abuse metrics but disqualify nobody")
    ap.add_argument("--no-mfg-gate", action="store_true",
                    help="report manufacturability but do not disqualify on it")
    ap.add_argument("--te-min-mm", type=float, default=None,
                    help=f"min hot-wire-cuttable thickness [mm] (default "
                         f"{CONFIG['te_min_mm']})")
    ap.add_argument("--gate-margin", type=float, default=None,
                    help=f"CLmax_gate must clear cl_target by this (default "
                         f"{CONFIG['gate_margin']})")
    ap.add_argument("--gate-ncrit", type=float, default=None,
                    help=f"ambient turbulence for the gate case (default "
                         f"{CONFIG['gate_n_crit']}; 9=avg tunnel, lower=windier)")
    ap.add_argument("--gate-xtr", type=float, default=None,
                    help=f"forced transition x/c for the gate case (default "
                         f"{CONFIG['gate_xtr']}; 0.05=fully tripped, 1.0=free)")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    if args.model_size:
        CONFIG["model_size"] = args.model_size       # now honored downstream too
    if args.gate_margin is not None:
        CONFIG["gate_margin"] = args.gate_margin
    if args.gate_ncrit is not None:
        CONFIG["gate_n_crit"] = args.gate_ncrit
    if args.gate_xtr is not None:
        CONFIG["gate_xtr"] = args.gate_xtr
    if args.te_min_mm is not None:
        CONFIG["te_min_mm"] = args.te_min_mm
    APPLY_ABUSE_GATE = not args.no_gate
    APPLY_MFG_GATE = not args.no_mfg_gate
    rng = np.random.default_rng(args.seed)

    print(f"MFG GATE ({'ON' if APPLY_MFG_GATE else 'OFF'}): hot-wire foam - section must "
          f"reach {CONFIG['te_min_mm']} mm thickness by "
          f"{100*(1-CONFIG['trunc_max_frac']):.0f}% chord "
          f"(truncating there costs <= {100*CONFIG['trunc_max_frac']:.0f}% of chord)")

    print(f"NeuralFoil model: {CONFIG['model_size']}")
    print(f"GATE  ({'ON' if APPLY_ABUSE_GATE else 'OFF'}): racing condition - "
          f"n_crit={CONFIG['gate_n_crit']}, xtr={CONFIG['gate_xtr']}, low-Re end.  "
          f"Need CLmax_gate >= cl_target + {CONFIG['gate_margin']}")
    print(f"ABUSE (reported, never gated): fully tripped LE, "
          f"xtr={CONFIG['abuse_xtr']}, n_crit={CONFIG['abuse_n_crit']}")
    print("Element profiles:")
    for pname, p in PROFILES.items():
        print(f"  {pname:14s} c={p['chord']*1000:.0f}mm  Re={p['re_list']}  "
              f"CL_target={p['cl_target']}  clmax_gate={p['clmax_gate']}  "
              f"t/c={p['tc']}  gate>={p['cl_target'] + CONFIG['gate_margin']:.2f}")

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
        p = PROFILES[pname]
        n_qual = int(df[col].notna().sum())
        print(f"\nTop 5 - {pname.upper()} (c={p['chord']*1000:.0f}mm, "
              f"{n_qual} qualifiers):")
        for _, r in top_n(df, col, 5).iterrows():
            ld = r[f"{pname}_LD_at_CL"]
            ld_s = "unreach" if np.isnan(ld) else f"{ld:6.1f}"
            print(f"  {r['name']:<22} CLmax={r[f'{pname}_CLmax']:.2f} "
                  f"L/D@CL{p['cl_target']}={ld_s} "
                  f"CLmax_gate={r[f'{pname}_CLmax_gate']:.2f} "
                  f"LD_ret_abuse={r[f'{pname}_LD_retention']:.2f} "
                  f"Cm={r[f'{pname}_Cm']:+.3f} score={r[col]:.3f}")


if __name__ == "__main__":
    main()
