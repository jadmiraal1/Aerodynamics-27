#!/usr/bin/env python3
"""
Shared config, geometry cache, and polar metrics for the front-wing screen.

Single source of truth for everything airfoil_screen.py, compare.py and
robustness.py must agree on. Before this module existed, each script kept its
own copy of the confidence thresholds and its own way of reducing a polar to
metrics, and they had already drifted apart (0.5 vs 0.6 valid-fraction; one
guarded against CLmax landing on the sweep boundary, the other did not).

Two rules this module exists to enforce:

1. RUN CONFIG IS DATA, NOT AN IMPORT.
   `from airfoil_screen import MODEL_SIZE` binds by value at import time, so a
   `--model-size` override in the screen never reached the downstream scripts:
   they would analyze a shortlist at a different fidelity than it was ranked at,
   silently. The screen now writes CONFIG to run_config.json next to
   results.xlsx, and compare/robustness load it back with load_run_config().

2. RANK AND RE-ANALYZE THE SAME GEOMETRY.
   The screen ranks Kulfan parameters. Downstream scripts used to rebuild
   airfoils from the Selig .dat files, which round-trips through a coordinate
   fit and yields a *slightly different section* than the one that was ranked.
   The screen now caches the exact Kulfan parameters to kulfan_params.json;
   load_kulfan() returns them, and NeuralFoil is called on those directly.
   The .dat files remain, but they are an XFLR5 export, not an input.
"""

import json
import os

import numpy as np

# ---------------------------------------------------------------------------
# Physical / operating point
# ---------------------------------------------------------------------------
NU = 1.46e-5                            # kinematic viscosity, sea level [m^2/s]
V_DESIGN = (16.0, 20.0, 27.0)           # lo / design / hi speeds [m/s] (lap sim)

ALPHA = np.arange(-2.0, 18.1, 0.5)      # AoA sweep, deg

# ---------------------------------------------------------------------------
# Element profiles. re_list is derived from chord x speeds (see below).
# ---------------------------------------------------------------------------
# cl_band is the OPERATING RANGE, and L/D is scored as the mean across it -
# not at a single point. Rationale: the wing's operating CL is set by INCIDENCE,
# which is a shim at the track, and the CL target itself inherits a wide
# uncertainty from the (still unfinished) chain lap sim -> ClA -> aero balance
# -> undertray share -> front wing load -> spanwise split. Scoring at one CL
# would let an undocumented constant pick the airfoil: sweeping main_center's
# target from 0.5 to 1.3 produced a DIFFERENT winning section at every step.
# Scoring across a band selects sections with a WIDE, FLAT drag bucket, whose
# ranking survives whatever the lap sim eventually says. Bands are +/-20% of
# nominal, reflecting the honest uncertainty in that chain.
PROFILES = {
    "main_center": dict(
        chord=0.275,
        cl_target=1.10,                 # NOMINAL: used for the gate and reporting.
                                        # See aero_targets.py - this is a judgment
                                        # inside a box the load budget narrows, NOT
                                        # a derived number. Bounded BELOW by what the
                                        # outboard can absorb, ABOVE by inlet wake.
        cl_band=(0.90, 1.35),           # operating range scored for L/D
        clmax_gate=1.30,                # still want stall margin above the band
        tc=(0.080, 0.160),              # same spar carries through
        weights={
            "LD_band":      0.40,       # wake thinness ACROSS the operating range
            "stall_gentle": 0.25,       # closest element to ground at the inlet
            "Cm_low":       0.20,       # front-loaded, not aft-loaded: less flow
                                        # turning into the nose. See README.
            "CL_usable":    0.05,
            "LD_max":       0.10,
        },
    ),
    "main_outboard": dict(
        chord=0.275,
        cl_target=1.50,
        cl_band=(1.25, 1.75),
        clmax_gate=1.40,                # hard gate on worst-case CLmax
        tc=(0.080, 0.160),              # spar needs thickness
        weights={
            "LD_band":      0.35,       # thin, low-loss wake -> undertray + balance
            "CL_usable":    0.15,
            "stall_gentle": 0.25,       # robust to ride height / pitch / yaw
            "Cm_low":       0.15,       # pitch sensitivity, structural twist
            "LD_max":       0.10,
        },
    ),
    "flap1": dict(
        chord=0.170,
        cl_target=1.80,
        cl_band=(1.55, 2.00),
        clmax_gate=1.60,
        tc=(0.055, 0.130),
        weights={
            "CLmax":        0.35,
            "CL_usable":    0.25,       # lift at 2 deg stall margin
            "stall_gentle": 0.20,       # tolerant to main-element upwash changes
            "LD_band":      0.10,
            "Cm_low":       0.10,       # relaxed: flap Cm reacts through the stack
        },
    ),
    "flap2": dict(
        chord=0.110,
        cl_target=1.75,                 # was 1.90: sat ~at the qualifier CLmax band,
                                        # so L/D was interpolated at the stall knee
                                        # (only 15% of qualifiers could reach it).
                                        # 1.75 -> ~50% reach, +0.10 median margin.
        cl_band=(1.50, 1.95),
        clmax_gate=1.60,
        tc=(0.045, 0.120),
        weights={
            "CLmax":        0.40,
            "CL_usable":    0.25,
            "stall_gentle": 0.20,
            "LD_band":      0.10,
            "Cm_low":       0.05,
        },
    ),
}
for _p in PROFILES.values():            # Re list from chord x (lo, design, hi)
    _p["re_list"] = [int(v * _p["chord"] / NU) for v in V_DESIGN]

TC_MIN_ALL = min(p["tc"][0] for p in PROFILES.values())
TC_MAX_ALL = max(p["tc"][1] for p in PROFILES.values())
ALL_RE = sorted({re for p in PROFILES.values() for re in p["re_list"]})
CL_TARGETS = {name: p["cl_target"] for name, p in PROFILES.items()}
CL_BANDS = {name: p["cl_band"] for name, p in PROFILES.items()}

# ---------------------------------------------------------------------------
# Run config. Mutable at runtime (CLI overrides), persisted to run_config.json,
# reloaded by the downstream scripts. Never import these values directly.
# ---------------------------------------------------------------------------
CONFIG = dict(
    # NeuralFoil ships 8 networks: xxsmall, xsmall, small, medium, large,
    # xlarge, xxlarge, xxxlarge. The library's own default is xlarge.
    #
    # This was "large" - BELOW the library default - with no justification. The
    # cost argument for staying small does not exist: a full 1,433-airfoil
    # screen takes 0.6 min at large and 3.3 min at xxxlarge. That is 2.7 minutes,
    # once per design cycle, to use the most accurate model available.
    #
    # And it matters. Re-scoring the shortlist across model sizes moves the
    # ranking (Spearman rho vs large: 0.61-0.78 for main_center, 0.96-0.99 for
    # the flaps), and main_center's winning FAMILY flips outright. The CLmax
    # differences between models are tiny (0.006-0.04) - the physics barely
    # moves. What moves is the ORDER, because the top candidates sit within
    # ~0.01 of each other in score, so model noise exceeds the gaps between
    # them. Use the best model and stop adding avoidable noise to a ranking
    # that is already this tightly packed.
    model_size="xxxlarge",
    conf_min=0.85,                      # min NeuralFoil analysis confidence
    min_valid_frac=0.6,                 # min fraction of alpha points passing conf
    camber_min=0.015,                   # symmetric sections pointless on a wing.
                                        # Tested: a symmetric section must run at
                                        # ~14 deg to make CL 1.1 at Re 300k, and
                                        # pays L/D 42 vs 78 for a cambered one -
                                        # a THICKER wake into the undertray, not a
                                        # cleaner one. See DESIGN_JUSTIFICATION.md.
    band_points=5,                      # CL samples across each element's cl_band

    # --- MANUFACTURABILITY: hot-wire cut foam cores ---
    # A knife-edge trailing edge cannot be cut in foam - the wire melts straight
    # through it and the edge crumbles. But TE thickness alone is the WRONG test:
    # a section can have an acceptable TE and still be a fragile feather at 90%
    # chord. So we ask the manufacturable question instead:
    #
    #   "At what chord station does this section reach a cuttable thickness,
    #    and can we afford to chop everything behind it?"
    #
    # Truncating at x_trunc gives a TE of exactly te_min_mm. The cost is the
    # chord you throw away - and with it, the aerodynamics you just ranked.
    # This matters enormously at flap chords: s1223 (the unconstrained flap2
    # winner) is 0.13 mm at the TE and still under 1 mm thick at 95% chord on a
    # 110 mm chord. Its thin cusped aft region is WHY it makes so much lift, and
    # it is exactly what cannot be manufactured.
    te_min_mm=1.5,                      # min cuttable thickness [mm]  <- SHOP LIMIT
    trunc_max_frac=0.03,                # max chord fraction we can afford to chop

    # --- THE GATE: the condition the car actually races in ---
    # Open-air autocross on a swept lot, with wind and cars ahead churning the
    # air, plus the dust/rubber a wing picks up over a 22 km endurance run.
    # n_crit 9 is Drela's "average wind tunnel", NOT free air; outdoor ambient
    # turbulence sits lower, and low-Re practice lands around 5-7. xtr 0.30
    # forces transition at 30% chord: a mildly dirty LE that still keeps some
    # laminar run - not a fully tripped one.
    #
    # Gating on the CLEAN case would be wrong in both directions. It is not
    # even the most favorable assumption: at these Re a long laminar run means
    # a fat laminar separation bubble, so forcing transition EARLIER often
    # RAISES CLmax. main_center qualifies 702 sections at n_crit 9 and 974 at
    # n_crit 6. A clean tunnel is a different flow, not a generous one.
    gate_n_crit=6.0,
    gate_xtr=0.30,
    gate_margin=0.00,                   # CLmax_gate must clear cl_target by this

    # --- THE ABUSE CASE: reported, never a disqualifier ---
    # Fully tripped LE (5% chord): rain, bugs, heavy rubber pickup. Worth
    # SEEING - a section that collapses here is fragile - but not worth
    # excluding on, since it is a scenario the car may never run in.
    abuse_xtr=0.05,
    abuse_n_crit=9.0,

    # --- robustness.py's three-case sweep ---
    track_n_crit=6.0,                   # realistic outdoor ambient turbulence
    clean_n_crit=9.0,                   # XFoil default / average tunnel
)

# Aft stations sampled for the truncation search. Starts at 70% because a
# section needing to be chopped that far forward is hopeless anyway.
AFT_X = np.linspace(0.70, 1.00, 31)

RUN_CONFIG_FILE = "run_config.json"
KULFAN_FILE = "kulfan_params.json"

# Transition environments. `clean` is the baseline the main screen ranks on;
# `abuse` is what the car actually sees once the LE picks up rubber.
CASES = {
    "clean": dict(n_crit=CONFIG["clean_n_crit"]),
    "track": dict(n_crit=CONFIG["track_n_crit"]),
    "abuse": dict(n_crit=CONFIG["clean_n_crit"],
                  xtr_upper=CONFIG["abuse_xtr"],
                  xtr_lower=CONFIG["abuse_xtr"]),
}


def gate_kwargs():
    """NeuralFoil kwargs for the GATE case - the condition the car races in."""
    return dict(n_crit=CONFIG["gate_n_crit"],
                xtr_upper=CONFIG["gate_xtr"],
                xtr_lower=CONFIG["gate_xtr"])


def abuse_kwargs():
    """NeuralFoil kwargs for the fully-tripped case. Reported, never gated on."""
    return dict(n_crit=CONFIG["abuse_n_crit"],
                xtr_upper=CONFIG["abuse_xtr"],
                xtr_lower=CONFIG["abuse_xtr"])


def save_run_config(outdir):
    with open(os.path.join(outdir, RUN_CONFIG_FILE), "w") as f:
        json.dump(CONFIG, f, indent=2)


def load_run_config(results_dir):
    """Load the config the results were actually produced with, into CONFIG.

    Downstream scripts MUST call this before analyzing, or they risk grading a
    shortlist at a fidelity it was never screened at.
    """
    path = os.path.join(results_dir, RUN_CONFIG_FILE)
    if not os.path.exists(path):
        print(f"  ! {path} not found - falling back to defaults. Results may not "
              f"match the run that produced {results_dir}.")
        return CONFIG
    with open(path) as f:
        CONFIG.update(json.load(f))
    return CONFIG


def save_kulfan(outdir, mapping):
    """mapping: name -> kulfan param dict (numpy arrays ok)."""
    ser = {n: {k: (v.tolist() if isinstance(v, np.ndarray) else v)
               for k, v in kp.items()}
           for n, kp in mapping.items()}
    with open(os.path.join(outdir, KULFAN_FILE), "w") as f:
        json.dump(ser, f)


def load_kulfan(results_dir):
    """Exact Kulfan params the screen ranked. Preferred over the .dat files."""
    path = os.path.join(results_dir, KULFAN_FILE)
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{path} not found. Re-run airfoil_screen.py to regenerate it - "
            f"re-fitting geometry from the .dat files analyzes a different "
            f"section than the one that was ranked.")
    with open(path) as f:
        raw = json.load(f)
    return {n: dict(upper_weights=np.array(d["upper_weights"]),
                    lower_weights=np.array(d["lower_weights"]),
                    leading_edge_weight=float(d["leading_edge_weight"]),
                    TE_thickness=float(d["TE_thickness"]))
            for n, d in raw.items()}


# ---------------------------------------------------------------------------
# Polar reduction - ONE implementation, used by every script.
# ---------------------------------------------------------------------------
def masked_polar(aero, alpha=ALPHA):
    """Confidence-masked (alpha, CL, CD, CM). None if the polar is unusable.

    Returns None when too few points clear CONFIG['conf_min'], or when CLmax
    lands on the sweep boundary (which means the true peak is outside the sweep
    and every downstream metric would be extrapolating).
    """
    conf = np.asarray(aero["analysis_confidence"], float)
    valid = conf >= CONFIG["conf_min"]
    if valid.mean() < CONFIG["min_valid_frac"]:
        return None
    a = np.asarray(alpha, float)[valid]
    cl = np.asarray(aero["CL"], float)[valid]
    cd = np.asarray(aero["CD"], float)[valid]
    cm = np.asarray(aero["CM"], float)[valid]

    i_max = int(np.argmax(cl))
    if i_max in (0, len(cl) - 1):       # CLmax on sweep boundary -> untrustworthy
        return None
    return a, cl, cd, cm, i_max


def ld_at(cl, cd, i_max, cl_target):
    """L/D at a target CL on the pre-stall branch.

    Returns np.nan (NOT 0.0) when the section cannot reach cl_target. A zero
    would be indistinguishable from "reached the target with catastrophic drag",
    and any later mean()/ratio over the column would quietly lie. Callers that
    want a score should map NaN -> 0 explicitly, at the point of scoring.
    """
    clmax = float(cl[i_max])
    if clmax < cl_target:
        return np.nan
    pre_cl, pre_cd = cl[:i_max + 1], cd[:i_max + 1]
    o = np.argsort(pre_cl)
    cd_t = float(np.interp(cl_target, pre_cl[o], pre_cd[o]))
    return cl_target / max(cd_t, 1e-6)


def normalized_metrics(sub, pname):
    """The six scored metrics, min-max normalized over the 5th-95th percentile
    of THIS element's qualifiers.

    Percentile clipping (rather than raw min-max) stops one outlier compressing
    everyone else into a narrow band. The cost is saturation at the top: several
    sections can tie at 1.0 on a metric.

    This is THE definition of a normalized metric, and both score() and
    weight_sensitivity.py call it. It used to live inside score(); the sensitivity
    study would have needed its own copy, and a copy is a bug waiting to happen -
    the study would eventually be measuring a scoring function the screen no
    longer used.
    """
    def norm(s, invert=False):
        lo, hi = s.quantile(0.05), s.quantile(0.95)
        x = ((s - lo) / max(hi - lo, 1e-9)).clip(0, 1)
        return 1 - x if invert else x

    return {
        "CLmax":        norm(sub[f"{pname}_CLmax"]),
        "CL_usable":    norm(sub[f"{pname}_CL_usable"]),
        "stall_gentle": norm(sub[f"{pname}_stall_gentle"]),
        "Cm_low":       norm(sub[f"{pname}_Cm"].abs(), invert=True),
        "LD_max":       norm(sub[f"{pname}_LD_max"]),
        "LD_band":      norm(sub[f"{pname}_LD_band"]),
    }


def truncation_station(aft_profile, chord_m):
    """Where must this section be chopped to give a cuttable trailing edge?

    aft_profile : local thickness / chord, sampled at AFT_X.
    chord_m     : this ELEMENT's chord [m] - the same section is manufacturable
                  at 275 mm and impossible at 110 mm, so this is per-element.

    Returns (x_trunc, chord_loss_frac, te_as_drawn_mm):
      x_trunc         furthest-aft station whose thickness >= te_min_mm.
                      Truncating here yields a TE of exactly te_min_mm.
      chord_loss_frac 1 - x_trunc: the chord you throw away to get there.
      te_as_drawn_mm  the section's OWN trailing-edge thickness, unmodified.

    If the section never reaches te_min_mm even at 70% chord, x_trunc is 0.70
    and chord_loss is 0.30 - which will fail any sane gate, as it should.
    """
    t_mm = np.asarray(aft_profile, float) * chord_m * 1000.0
    ok = np.nonzero(t_mm >= CONFIG["te_min_mm"])[0]
    if len(ok) == 0:
        return float(AFT_X[0]), float(1.0 - AFT_X[0]), float(t_mm[-1])
    x_trunc = float(AFT_X[ok[-1]])          # furthest-aft cuttable station
    return x_trunc, float(1.0 - x_trunc), float(t_mm[-1])


def ld_band_mean(cl, cd, i_max, band, n_points=None):
    """Mean L/D across the operating CL band. THIS IS WHAT GETS SCORED.

    A CL inside the band that the section cannot reach scores ZERO, not NaN -
    and that is deliberate, unlike ld_at(). Here "unreachable" is a real,
    quantified deficiency: the wing is expected to operate at that CL, and this
    section cannot deliver it. Averaging zeros in penalizes exactly in
    proportion to how much of the operating range the section forfeits. A
    section that covers the whole band beats one that is brilliant over half of
    it and absent over the rest - which is the entire point of band scoring.
    """
    n = n_points or CONFIG["band_points"]
    vals = [ld_at(cl, cd, i_max, t) for t in np.linspace(band[0], band[1], n)]
    return float(np.mean([0.0 if np.isnan(v) else v for v in vals]))


def ld_band_flat(cl, cd, i_max, band, n_points=None):
    """Bucket flatness across the band: min/max of L/D, in [0, 1].

    1.0 = perfectly flat drag bucket across the whole operating range (the
    section does not care where you trim it). Low = the section has a sharp
    optimum and falls off either side, so its ranking depends on a CL target
    you do not actually know yet. Reported, not scored - it is a confidence
    measure on the ranking, not a performance metric.
    """
    n = n_points or CONFIG["band_points"]
    vals = [ld_at(cl, cd, i_max, t) for t in np.linspace(band[0], band[1], n)]
    vals = [0.0 if np.isnan(v) else v for v in vals]
    hi = max(vals)
    return float(min(vals) / hi) if hi > 0 else 0.0


def polar_metrics(aero, alpha=ALPHA, cl_targets=None):
    """Full metric set from one polar (one Re, one transition case)."""
    m = masked_polar(aero, alpha)
    if m is None:
        return None
    a, cl, cd, cm, i_max = m
    clmax, a_stall = float(cl[i_max]), float(a[i_max])

    a_use = a_stall - 2.0               # usable point: 2 deg below stall
    cl_use = float(np.interp(a_use, a, cl))
    cm_use = float(np.interp(a_use, a, cm))

    a_post = min(a_stall + 3.0, a[-1])  # CL retention 3 deg past stall (0..1)
    gentle = max(0.0, min(1.0, float(np.interp(a_post, a, cl)) / clmax))

    cl_pre, cd_pre = cl[:i_max + 1], cd[:i_max + 1]
    targets = CL_TARGETS if cl_targets is None else cl_targets
    return dict(
        CLmax=clmax, alpha_stall=a_stall, CL_usable=cl_use,
        stall_gentle=gentle, Cm_use=cm_use,
        LD_max=(float(np.max(cl_pre[1:] / np.maximum(cd_pre[1:], 1e-6)))
                if i_max > 1 else 0.0),
        # single-point L/D at the nominal target: REPORTED
        LD_at={name: ld_at(cl, cd, i_max, t) for name, t in targets.items()},
        # mean L/D across the operating band: SCORED
        LD_band={name: ld_band_mean(cl, cd, i_max, CL_BANDS[name])
                 for name in CL_BANDS},
        LD_flat={name: ld_band_flat(cl, cd, i_max, CL_BANDS[name])
                 for name in CL_BANDS},
    )
