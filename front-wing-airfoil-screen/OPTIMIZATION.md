# Rigging Optimization Study — `optimize.py`

2-D multi-element rigging optimization in ground effect, driving SU2.

## The objective

```
minimize  Cd     subject to   CL_down >= 2.03
```

**Minimum drag AT the required downforce — not maximum downforce.** `aero_targets.py` already told us how much the front wing must make (187 N at 20 m/s, the front axle's share of ClA = 4.0). Downforce past that target does not help the car — it *unbalances* it, and you pay for it in drag. So the optimizer is asked to hit the number and spend as little drag as possible doing so.

## Read this before you run 30 CFD cases

**With 20–50 runs over 7 variables you cannot find a global optimum, and you cannot fit a surrogate worth trusting.** An 8-dimensional space needs hundreds of points before a Gaussian process means anything. Anyone claiming a global optimum from 30 runs is fooling themselves, and a design judge will find the seam.

What 30 runs *can* buy you is **which variables matter** — and then a real optimization over the two or three that do. That is what this does.

## The three phases

### Phase 1 — Morris screening (elementary effects)

`r` trajectories × (k+1) runs. Default r=3, k=7 → **24 runs**. The cheapest statistically defensible way to rank 7 variables.

| Output | Meaning |
|---|---|
| **μ\*** | how much this variable moves the objective |
| **σ** | how much its effect *depends on the others* — i.e. interaction strength |

A variable with high μ\* and low σ can be tuned on its own. High σ means it interacts and must be tuned jointly. Both are things you want to know *before* spending your remaining budget.

### Phase 2 — Local refinement

Nelder-Mead over **only** the top-N variables from Phase 1, holding the rest at their screened-best. ~10–20 runs.

### Phase 3 — Ride-height robustness

**Ride height is an operating condition, not a design variable** — you cannot bolt a ride height to the car. It pitches and heaves. A rigging that is optimal at 50 mm and falls apart at 30 mm is not optimal, it is *lucky*.

The finalists are re-run across the ride-height band (30 / 50 / 70 mm) and ranked on **worst-case** objective, not mean.

## Usage

```
python optimize.py --baseline                    # ONE case. Do this FIRST.
python optimize.py --solver dummy --phase all    # test the harness, no CFD
python optimize.py --phase screen                # Morris, 24 runs
python optimize.py --phase refine --n-active 3   # local search
python optimize.py --phase robust                # ride-height sweep
```

Outputs land in `<results>/opt_study/`: `study.csv` (every evaluation), `screening.csv` (μ\*/σ), `robustness.csv`, and `cases/` (one SU2 case dir each).

## Run the baseline first. Seriously.

```
python optimize.py --baseline
```

**Never spend a DOE budget on an unvalidated solver setup.** Check three things:

1. It converges.
2. Forces are physically sane.
3. **CL is NEGATIVE.**

That last one matters: `mesh_section.py` mirrors the stack into car orientation, so the wing pushes *down*. SU2 reports CL positive-up, so a working front wing gives a **negative CL**. The code optimizes `cl_down = -CL`. Get this backwards and you will run a beautiful, expensive study that carefully maximizes *lift*.

## Everything is checkpointed

Every evaluation is written to `study.csv` immediately, keyed on the rigging + ride height. Kill it, restart it, extend it when HPC access appears — **it never repeats a run it has already paid for.** Re-running a completed phase costs zero.

This also means widening a bound and re-running is cheap: only the new points cost anything.

## Guards that save runs (and save you)

**Geometric feasibility, checked before meshing.** Many rigging combinations are physically impossible (elements intersecting) or unmeshable (slot so tight the BL extrusions collide). These are rejected in milliseconds instead of costing a mesh + a CFD run. Verified: slots at 0.4% chord → rejected (0.7 mm gap); ride height 2 mm → rejected.

**Bound-hitting warning.** If the best design sits *on* one of its bounds, the optimizer says so, loudly. **This is the easiest way to ship a wrong answer from a DOE**: the optimizer cannot tell you "I would have gone further" — it just stops at the edge, and the result *looks* like a converged optimum. It isn't. It is the edge of where you let it look. Widen `DESIGN_VARS` and re-run (the study resumes).

## The weakest number in the study

```python
CL_TARGET = 2.03
```

`aero_targets.py` gives the front wing a **3-D** CL of ~1.74 across its span, which — with a 30% main-plane-only centre — implies the outboard stack must make CL ≈ 2.03.

**That is a 3-D wing CL. This is a 2-D section simulation.** They are not the same quantity: 2-D omits tip losses and the endplate/vortex system that dominates a real FSAE front wing. 2.03 is a placeholder with the right order of magnitude, **not a derived target**.

Fix it with a 3-D correction factor as soon as you have one, and re-run — **the optimal rigging depends on the loading you ask for.** Ask for the wrong loading and you get the right rigging for the wrong wing.

## What is *not* modeled

- **3-D everything.** No endplates, no tip vortices, no outwash. On a real FSAE front wing these are first-order, not a correction. A 2-D rigging optimum is a *starting point* for 3-D CFD, not an answer.
- **The centre section.** This meshes the outboard 3-element stack only. The centre (main plane alone, feeding the undertray) is a different section and a different problem.
- **Wheel and chassis interaction.** The front wing's whole downstream job — feeding the undertray, managing wheel wake — is invisible to a 2-D section.

## Validation status — read this honestly

**The harness is fully tested. The SU2 setup is not.**

I could not run SU2 (or even gmsh) in the environment this was written in. So:

- ✅ **Harness:** validated end to end against `--solver dummy`, an analytic stand-in with a real optimum, realistic variable interactions, and a stall-like penalty. All three phases run, resume works (24 cached / 0 re-run), feasibility rejection fires, bound-warnings fire.
- ⚠️ **SU2 config:** written from documentation, not from a run that converged. SU2 config keys drift between major versions. **`--baseline` is not optional.**

The dummy solver exists precisely so that every bug in the *search logic* is found in 2 seconds rather than 6 hours into a DOE. Use it whenever you change the harness.

## Suggested first session

```
python optimize.py --solver dummy --phase all     # 30 s. Proves the harness.
python optimize.py --baseline                     # 1 CFD run. Proves SU2.
python optimize.py --phase screen                 # 24 runs. What matters?
python optimize.py --phase refine --n-active 3    # ~18 runs.
python optimize.py --phase robust                 # 9 runs. Does it survive pitch?
```

~52 CFD runs. Stop after screening if the answer is already obvious — the μ\*/σ table often tells you most of what you need.
