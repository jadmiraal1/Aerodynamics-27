#!/usr/bin/env python3
"""
Rigging optimization study: 2-D multi-element front wing in ground effect.

OBJECTIVE
---------
    minimize  Cd     subject to   Cd_down >= CL_TARGET

Minimum drag AT the required downforce - not maximum downforce. The load budget
(aero_targets.py) already told us how much downforce the front wing must make:
187 N at 20 m/s, i.e. the front axle's share of a ClA of 4.0. Downforce beyond
that target does not help the car - it unbalances it, and it is paid for in
drag. So the optimizer is asked to HIT the number and spend as little drag as
possible doing it.

Constraint handled by penalty, so infeasible designs still carry gradient
information and the search can walk back into the feasible region rather than
falling off a cliff.

BE HONEST ABOUT THE BUDGET
-------------------------
With ~20-50 CFD runs over 7 rigging variables, you CANNOT find a global optimum
and you cannot fit a trustworthy surrogate. 8-dimensional space needs hundreds
of points before a Gaussian process means anything. Anyone claiming a global
optimum from 30 runs is fooling themselves, and a design judge will know.

What you CAN do with 30 runs is find out WHICH VARIABLES MATTER, and then
optimize over the two or three that do. That is what this does:

  PHASE 1 - SCREENING (Morris elementary effects).
      r trajectories x (k+1) runs = the cheapest statistically defensible way
      to rank 7 variables. Gives mu* (how much each variable moves the answer)
      and sigma (how much its effect depends on the others, i.e. interactions).
      Default r=3 -> 24 runs.

  PHASE 2 - REFINEMENT.
      Local search (Nelder-Mead) over only the variables Phase 1 says matter,
      holding the rest at their screened-best values. ~10-20 runs.

  PHASE 3 - RIDE-HEIGHT ROBUSTNESS.
      The car pitches and heaves; a rigging that is optimal at one ride height
      and falls apart 10 mm away is not optimal, it is lucky. The finalists are
      re-run across the ride-height band. Ride height is an OPERATING
      CONDITION, not a design variable - you cannot bolt a ride height to the
      car.

EVERY EVALUATION IS CHECKPOINTED to study.csv. Kill it, restart it, extend it
when HPC access appears - it never repeats a run it has already paid for.

Usage:
  python optimize.py --baseline                 # ONE case. Do this FIRST.
  python optimize.py --solver dummy --phase all # test the harness (no CFD)
  python optimize.py --phase screen             # Morris screening
  python optimize.py --phase refine             # local search on what matters
  python optimize.py --phase robust             # ride-height sweep of finalists
"""

import argparse
import itertools
import json
import os
import time

import numpy as np
import pandas as pd

from common import PROFILES, V_DESIGN, load_run_config
from solvers import get_solver
import mesh_section as MS

# ---------------------------------------------------------------------------
# DESIGN SPACE
#
# Bounds are engineering judgment, not derived. They bracket what is riggable
# on the car and meshable by mesh_section.py. Widen them and you will spend
# CFD runs on geometry the mesher rejects; narrow them and you may bracket out
# the optimum.
# ---------------------------------------------------------------------------
DESIGN_VARS = {
    "alpha_main": (-2.0, 8.0),      # main incidence [deg]
    "delta1":     (10.0, 32.0),     # flap1 deflection vs main chord [deg]
    "delta2":     (12.0, 40.0),     # flap2 deflection vs flap1 chord [deg]
    "gap1":       (0.012, 0.040),   # slot gap, fraction of parent chord
    "gap2":       (0.012, 0.040),
    "overlap1":   (0.000, 0.060),   # flap LE upstream of parent TE
    "overlap2":   (0.000, 0.060),
}

# Ride height is an OPERATING CONDITION, not a design variable.
RIDE_HEIGHTS = [0.030, 0.050, 0.070]        # [m]
RIDE_NOMINAL = 0.050

# Target downforce coefficient (2-D, per unit span).
#
# !! THIS IS THE WEAKEST NUMBER IN THE STUDY. !!
# aero_targets.py gives the front wing a 3-D CL of ~1.74 over its whole span,
# which with a 30% main-plane-only center implies the OUTBOARD stack must make
# CL ~2.03. That is a 3-D wing CL. This is a 2-D section simulation. They are
# not the same quantity: 2-D omits tip losses and the endplate/vortex system
# that dominates a real FSAE front wing.
#
# 2.03 is therefore a placeholder with the right order of magnitude, not a
# derived target. Fix it with a 3-D CFD correction factor as soon as you have
# one, and re-run - the optimum rigging depends on the loading you ask for.
CL_TARGET = 2.03
PENALTY = 5.0                        # drag counts per unit of missed downforce


def feasible(rig):
    """Cheap geometric rejection - BEFORE spending a CFD run.

    Many rigging combinations produce geometry that is physically impossible
    (elements intersecting) or unmeshable (a slot so tight the BL extrusions
    collide). Meshing then CFD-ing them wastes the scarcest thing in this study.
    Returns (ok, reason).
    """
    secs = _sections()
    try:
        placed, _ = MS.rig(secs, dict(MS.RIG_DEFAULTS, **rig))
    except Exception as e:
        return False, f"rig failed: {e}"

    g1 = MS.min_gap(placed["main"], placed["flap1"])
    g2 = MS.min_gap(placed["flap1"], placed["flap2"])
    if min(g1, g2) < 2.0e-3:
        return False, f"slot gap {min(g1,g2)*1000:.1f} mm < 2.0 mm (unmeshable)"

    ground = min(p[:, 1].min() for p in placed.values())
    if ground < 5.0e-3:
        return False, f"ground clearance {ground*1000:.1f} mm < 5 mm"
    return True, ""


_SEC_CACHE = {}


def _sections():
    """Blunted rank-1 sections. Cached - identical for every design."""
    if "s" not in _SEC_CACHE:
        raise RuntimeError("call load_sections() first")
    return _SEC_CACHE["s"]


def load_sections(results_dir):
    load_run_config(results_dir)
    names = MS.rank1_sections(results_dir)
    secs = {}
    for (elem, profile, npts), name in zip(MS.ELEMENTS, names):
        af = MS.load_section(name, results_dir, profile_key=profile)
        secs[elem] = MS.unit_coords(af, npts)
    _SEC_CACHE["s"] = secs
    _SEC_CACHE["names"] = names
    return names


def objective(res):
    """min Cd s.t. cl_down >= CL_TARGET, as a single penalized number."""
    if not res["converged"]:
        return np.inf
    shortfall = max(0.0, CL_TARGET - res["cl_down"])
    return res["cd"] + PENALTY * shortfall ** 2


# ---------------------------------------------------------------------------
# Evaluation with checkpointing
# ---------------------------------------------------------------------------
class Study:
    def __init__(self, outdir, solver, results_dir, v_inf):
        self.outdir = outdir
        self.csv = os.path.join(outdir, "study.csv")
        self.solver = solver
        self.results_dir = results_dir
        self.v_inf = v_inf
        os.makedirs(outdir, exist_ok=True)
        self.df = (pd.read_csv(self.csv) if os.path.exists(self.csv)
                   else pd.DataFrame())
        if len(self.df):
            print(f"  resuming: {len(self.df)} evaluations already on disk")

    def _key(self, rig, h):
        return json.dumps({**{k: round(float(rig[k]), 6) for k in DESIGN_VARS},
                           "ride_height": round(float(h), 6)}, sort_keys=True)

    def evaluate(self, rig, h=RIDE_NOMINAL, tag=""):
        key = self._key(rig, h)
        if len(self.df) and (self.df["key"] == key).any():
            row = self.df[self.df["key"] == key].iloc[-1]
            return dict(row), True                       # cached: costs nothing

        full = dict(MS.RIG_DEFAULTS, **rig, ride_height=h)
        ok, why = feasible(full)
        if not ok:
            res = dict(cl=np.nan, cd=np.nan, cl_down=np.nan, converged=False,
                       error=why)
        elif self.solver.name == "dummy":
            # The dummy solver has no mesh. Building one would defeat the whole
            # point of it - it exists so the HARNESS can be exercised in seconds
            # without gmsh, without a solver, and without spending CFD runs on
            # bugs in the search logic. It still goes through feasible() above,
            # so infeasible-design rejection is exercised for real.
            t0 = time.time()
            res = self.solver.run(None, None, self.v_inf, self.ref_chord, rig=full)
            res["wall_s"] = time.time() - t0
        else:
            case = os.path.join(self.outdir, "cases",
                                f"{len(self.df):04d}{('_' + tag) if tag else ''}")
            mesh = self._mesh(full, case)
            if mesh is None:
                res = dict(cl=np.nan, cd=np.nan, cl_down=np.nan,
                           converged=False, error="mesh failed")
            else:
                t0 = time.time()
                res = self.solver.run(mesh, case, self.v_inf, self.ref_chord)
                res["wall_s"] = time.time() - t0

        row = {**{k: full[k] for k in DESIGN_VARS}, "ride_height": h,
               "key": key, "tag": tag,
               **{k: res.get(k) for k in
                  ("cl", "cd", "cl_down", "converged", "iters", "error", "wall_s")}}
        row["objective"] = objective(res)
        self.df = pd.concat([self.df, pd.DataFrame([row])], ignore_index=True)
        self.df.to_csv(self.csv, index=False)            # checkpoint EVERY run
        return row, False

    ref_chord = 0.275

    def _mesh(self, full, case_dir):
        """Build the mesh for one rigging. Returns the .su2 path, or None."""
        os.makedirs(case_dir, exist_ok=True)
        try:
            placed, info = MS.rig(_sections(), full)
            gaps = dict(main_flap1=MS.min_gap(placed["main"], placed["flap1"]),
                        flap1_flap2=MS.min_gap(placed["flap1"], placed["flap2"]))
            spec = MS.bl_spec(self.v_inf, 1.0, min(gaps.values()))
            c = PROFILES["main_outboard"]["chord"]
            bx1 = max(p[:, 0].max() for p in placed.values())
            by1 = max(p[:, 1].max() for p in placed.values())
            dom = dict(x0=-8 * c, x1=16 * c, y0=0.0, y1=10 * c,
                       far_size=0.15, near_size=2.0e-3, near_dist=0.10,
                       wake_size=6.0e-3, wake_x0=-0.05, wake_x1=bx1 + 4 * c,
                       wake_y1=by1 + 0.05,
                       gnd_size=2.5e-3, gnd_x0=-0.3, gnd_x1=bx1 + 3 * c,
                       gnd_h=0.02)
            base = os.path.join(case_dir, "mesh")
            MS.build_mesh(placed, _sections(), spec, dom, base, ["su2"])
            return base + ".su2"
        except Exception as e:
            print(f"    mesh failed: {e}")
            return None


# ---------------------------------------------------------------------------
# PHASE 1 - Morris screening
# ---------------------------------------------------------------------------
def morris(study, r=3, levels=4, seed=0):
    """Elementary effects. Cost = r*(k+1) runs. Returns mu*, sigma per variable.

    mu*  = mean |effect| -> how much this variable moves the objective.
    sigma= spread of the effect -> how much it DEPENDS ON the others.
           High sigma means the variable interacts and cannot be tuned alone.
    """
    keys = list(DESIGN_VARS)
    k = len(keys)
    rng = np.random.default_rng(seed)
    delta = levels / (2.0 * (levels - 1))
    lo = np.array([DESIGN_VARS[x][0] for x in keys])
    hi = np.array([DESIGN_VARS[x][1] for x in keys])

    effects = {x: [] for x in keys}
    n = r * (k + 1)
    print(f"\nPHASE 1 - Morris screening: {r} trajectories x {k+1} = {n} runs")

    def to_rig(u):
        v = lo + u * (hi - lo)
        return {x: float(v[i]) for i, x in enumerate(keys)}

    done = 0
    for t in range(r):
        base = rng.integers(0, levels // 2, size=k) / (levels - 1.0)
        order = rng.permutation(k)
        u = base.copy()
        row, cached = study.evaluate(to_rig(u), tag=f"morris_t{t}_0")
        done += 1
        prev = row["objective"]
        print(f"  [{done}/{n}] base       obj={prev:.5f}"
              + ("  (cached)" if cached else ""))
        for j, i in enumerate(order):
            u2 = u.copy()
            u2[i] = min(1.0, u2[i] + delta) if u2[i] + delta <= 1.0 else u2[i] - delta
            row, cached = study.evaluate(to_rig(u2), tag=f"morris_t{t}_{keys[i]}")
            done += 1
            cur = row["objective"]
            if np.isfinite(cur) and np.isfinite(prev):
                effects[keys[i]].append((cur - prev) / delta)
            print(f"  [{done}/{n}] {keys[i]:<11s} obj={cur:.5f}"
                  + ("  (cached)" if cached else ""))
            u, prev = u2, cur

    out = []
    for x in keys:
        e = np.array(effects[x], float)
        out.append(dict(var=x,
                        mu_star=float(np.mean(np.abs(e))) if len(e) else np.nan,
                        sigma=float(np.std(e)) if len(e) else np.nan,
                        n=len(e)))
    return pd.DataFrame(out).sort_values("mu_star", ascending=False)


# ---------------------------------------------------------------------------
# PHASE 2 - local refinement on the variables that matter
# ---------------------------------------------------------------------------
def refine(study, active, x0, budget=20):
    """Nelder-Mead over `active` vars only; the rest stay at x0."""
    from scipy.optimize import minimize

    lo = np.array([DESIGN_VARS[x][0] for x in active])
    hi = np.array([DESIGN_VARS[x][1] for x in active])
    calls = {"n": 0}

    def f(u):
        if calls["n"] >= budget:
            return np.inf
        v = np.clip(u, 0, 1) * (hi - lo) + lo
        rig = dict(x0, **{a: float(v[i]) for i, a in enumerate(active)})
        row, cached = study.evaluate(rig, tag="refine")
        if not cached:
            calls["n"] += 1
            print(f"  [{calls['n']}/{budget}] "
                  + "  ".join(f"{a}={rig[a]:.4g}" for a in active)
                  + f"   obj={row['objective']:.5f}")
        return row["objective"] if np.isfinite(row["objective"]) else 1e3

    u0 = np.array([(x0[a] - DESIGN_VARS[a][0])
                   / (DESIGN_VARS[a][1] - DESIGN_VARS[a][0]) for a in active])
    print(f"\nPHASE 2 - refining {active} ({budget} runs max)")
    minimize(f, u0, method="Nelder-Mead",
             options=dict(maxfev=budget, xatol=0.02, fatol=1e-4))
    return study


# ---------------------------------------------------------------------------
# PHASE 3 - ride-height robustness
# ---------------------------------------------------------------------------
def robustness(study, finalists):
    print(f"\nPHASE 3 - ride-height robustness ({len(finalists)} designs "
          f"x {len(RIDE_HEIGHTS)} heights)")
    rows = []
    for i, rig in enumerate(finalists):
        for h in RIDE_HEIGHTS:
            row, _ = study.evaluate(rig, h=h, tag=f"robust{i}")
            rows.append(dict(design=i, ride_height=h, cl_down=row["cl_down"],
                             cd=row["cd"], objective=row["objective"]))
            print(f"  design {i}  h={h*1000:.0f} mm  "
                  f"CL_down={row['cl_down']:.3f}  Cd={row['cd']:.4f}")
    d = pd.DataFrame(rows)
    # A rigging that only works at one ride height is lucky, not optimal.
    summ = d.groupby("design").agg(
        cl_mean=("cl_down", "mean"), cl_spread=("cl_down", lambda s: s.max()-s.min()),
        cd_mean=("cd", "mean"), obj_worst=("objective", "max")).reset_index()
    return d, summ.sort_values("obj_worst")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="screen_results")
    ap.add_argument("--out", default=None, help="default <results>/opt_study")
    ap.add_argument("--solver", default="su2", choices=["su2", "dummy"])
    ap.add_argument("--phase", default="all",
                    choices=["all", "screen", "refine", "robust"])
    ap.add_argument("--baseline", action="store_true",
                    help="run ONE case at the default rigging and stop")
    ap.add_argument("--trajectories", type=int, default=3)
    ap.add_argument("--refine-budget", type=int, default=20)
    ap.add_argument("--n-active", type=int, default=3,
                    help="how many top variables to refine over")
    ap.add_argument("--v", type=float, default=V_DESIGN[1],
                    help=f"freestream [m/s], default design speed {V_DESIGN[1]}")
    ap.add_argument("--su2-exe", default="SU2_CFD")
    ap.add_argument("--iters", type=int, default=3000)
    ap.add_argument("--transition", action="store_true",
                    help="Langtry-Menter transition model (slower, consistent "
                         "with the screen's physics)")
    ap.add_argument("--nproc", type=int, default=1)
    args = ap.parse_args()

    outdir = args.out or os.path.join(args.results, "opt_study")
    names = load_sections(args.results)
    print(f"Sections: {dict(zip([e for e,_,_ in MS.ELEMENTS], names))}")
    print(f"Objective: min Cd  s.t.  CL_down >= {CL_TARGET}  (penalty {PENALTY})")
    print(f"Freestream: {args.v} m/s   solver: {args.solver}")

    solver = get_solver(args.solver, exe=args.su2_exe, iters=args.iters,
                        transition=args.transition, nproc=args.nproc) \
        if args.solver == "su2" else get_solver("dummy")
    if args.solver == "su2" and not solver.available():
        raise SystemExit(f"'{args.su2_exe}' not on PATH. Install SU2 or use "
                         f"--solver dummy to test the harness.")

    study = Study(outdir, solver, args.results, args.v)

    if args.baseline:
        print("\nBASELINE - one case at the default rigging.")
        print("Check it converges and that CL is NEGATIVE (front wing pushes "
              "down) before spending a DOE on it.")
        row, _ = study.evaluate(dict(MS.RIG_DEFAULTS), tag="baseline")
        print(f"\n  CL (SU2 sign) = {row['cl']:.4f}   <- must be NEGATIVE")
        print(f"  CL_down       = {row['cl_down']:.4f}   (target {CL_TARGET})")
        print(f"  Cd            = {row['cd']:.5f}")
        print(f"  objective     = {row['objective']:.5f}")
        print(f"  case: {outdir}/cases/")
        return

    screen_csv = os.path.join(outdir, "screening.csv")
    if args.phase in ("all", "screen"):
        eff = morris(study, r=args.trajectories)
        eff.to_csv(screen_csv, index=False)
        print("\n  MORRIS ELEMENTARY EFFECTS")
        print("  mu*   = how much the variable moves the objective")
        print("  sigma = how much its effect depends on the others (interaction)")
        print(eff.to_string(index=False))

    if args.phase in ("all", "refine"):
        if not os.path.exists(screen_csv):
            raise SystemExit("run --phase screen first")
        eff = pd.read_csv(screen_csv)
        active = eff.head(args.n_active)["var"].tolist()
        best = study.df.loc[study.df["objective"].idxmin()]
        x0 = {k: float(best[k]) for k in DESIGN_VARS}
        print(f"\n  refining over: {active}")
        print(f"  holding fixed: {[k for k in DESIGN_VARS if k not in active]}")
        refine(study, active, x0, budget=args.refine_budget)

    if args.phase in ("all", "robust"):
        d = study.df[study.df["converged"] == True]  # noqa: E712
        d = d[d["ride_height"] == RIDE_NOMINAL].nsmallest(3, "objective")
        finalists = [{k: float(r[k]) for k in DESIGN_VARS}
                     for _, r in d.iterrows()]
        det, summ = robustness(study, finalists)
        det.to_csv(os.path.join(outdir, "robustness.csv"), index=False)
        print("\n  RIDE-HEIGHT ROBUSTNESS (worst-case objective wins)")
        print(summ.to_string(index=False))

    ok = study.df[study.df["converged"] == True]  # noqa: E712
    if len(ok):
        b = ok.loc[ok["objective"].idxmin()]
        print("\n" + "=" * 68)
        print("BEST DESIGN")
        print("=" * 68)
        on_bound = []
        for k, (lo, hi) in DESIGN_VARS.items():
            v = float(b[k])
            edge = ("  <-- LOWER BOUND" if v - lo < 0.02 * (hi - lo) else
                    "  <-- UPPER BOUND" if hi - v < 0.02 * (hi - lo) else "")
            if edge:
                on_bound.append(k)
            print(f"  {k:<12s} {v:8.4g}   [{lo}, {hi}]{edge}")
        print(f"  -> CL_down {b['cl_down']:.3f} (target {CL_TARGET}), "
              f"Cd {b['cd']:.5f}")

        # A design pressed against its box is telling you the box is wrong. The
        # optimizer cannot say "I would have gone further" - it just stops, and
        # the result LOOKS like an optimum. This is the single easiest way to
        # ship a wrong answer from a DOE, so say it loudly.
        if on_bound:
            print(f"\n  !! {len(on_bound)} variable(s) sit ON their bounds: "
                  f"{', '.join(on_bound)}")
            print("     The optimum is probably OUTSIDE the design box. These are")
            print("     not converged values - they are the edge of where you let")
            print("     the optimizer look. Widen DESIGN_VARS for these and re-run")
            print("     (the study resumes; nothing already paid for is repeated).")

        n_fail = int((study.df["converged"] != True).sum())  # noqa: E712
        print(f"\n  {len(study.df)} evaluations, {n_fail} infeasible/failed")
        print(f"  all runs: {study.csv}")


if __name__ == "__main__":
    main()
