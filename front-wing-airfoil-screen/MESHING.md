# 2D CFD Meshing Stage — `mesh_section.py`

Rigs shortlist sections into the outboard three-element stack and produces a
solver-ready 2D RANS mesh in ground effect. Solver-agnostic on purpose: this
year's compute/licensing is not locked in, so the mesh must not lock it in
either.

## Quick start

```
pip install -r requirements.txt                 # gmsh is included
python run_pipeline.py                          # meshing is stage 5, runs by default

python mesh_section.py                          # or run the stage alone:
python mesh_section.py --airfoils MAIN F1 F2    #   explicit sections
python mesh_section.py --ride-height 0.035 --delta1 24 --gap1 0.02
```

**Linux only:** the gmsh wheel links against system OpenGL libs pip cannot install. If `import gmsh` raises `libGLU.so.1: cannot open shared object file`, run `sudo apt-get install libglu1-mesa libxcursor1 libxinerama1 libxft2`. Windows and macOS wheels bundle these.

## The trailing edge is blunted before meshing — and it must be

The screen does **not** rank the airfoil as drawn. These wings are hot-wire cut from foam, so a knife-edge TE cannot be manufactured; `airfoil_screen.py` sets `TE_thickness = te_min_mm / chord` **per element** before running any polar, and ranks the *blunted* section.

That blunt is per-element, because it is an absolute thickness: 1.5 mm is 0.55% of the 275 mm main chord but **1.36%** of the 110 mm flap2 chord. A single per-name geometry file cannot carry both, which is why `kulfan_params.json` stores the **as-drawn** section and every consumer re-applies the blunt itself. `mesh_section.py` reads `te_min_mm` from the run's own `run_config.json` and does exactly that:

```
goe63_v16 [main_outboard]: TE 1.93 mm -> 1.93 mm   (already thick enough)
s1223_v19 [flap1]        : TE 0.21 mm -> 1.50 mm
s1223_v19 [flap2]        : TE 0.13 mm -> 1.50 mm
```

Mesh the as-drawn geometry instead and you solve a knife-edge section that was never screened, cannot be built, and whose aft loading is not the one that won.

Note `TE_COLLAPSE` (3e-4, unit chord) collapses *sharper* bases to a point to avoid micro-elements. A 1.5 mm blunt at flap2 is 1.36e-2 unit chord — far above it, so the blunt survives and the TE base is meshed as a real edge.

Outputs land in **`<results>/mesh_out/`** — i.e. `screen_results/mesh_out/` by default, *inside* the results directory, next to `results.xlsx` and `run_config.json`. The mesh belongs with the run that produced it: a cwd-relative `mesh_out/` had no link back to the config it was built from, so a re-screen would silently leave a stale mesh sitting beside fresh results. `--out DIR` overrides.

Each run writes: the mesh in five formats + `*_report.json` (counts, quality, y+ basis, measured gaps, BL spec) + `*_preview.png`. **Always look at the preview before solving.**

## Mesh settings

Defaults, for the rank-1 stack at 27 m/s:

| Setting | Value | Note |
|---|---|---|
| `--yplus` | **1.0** | transition-SST ready. The entire screen rests on laminar/transition physics — wall functions would throw that away. |
| `--v-ref` | `V_DESIGN[-1]` = 27 m/s | **imported from `common.py`, not hardcoded.** Highest Re → thinnest first cell → conservative for every lower speed. |
| first cell h₁ | 11.8 µm | from a turbulent flat-plate estimate at the main-element Re (508k) |
| growth ratio | 1.2 | |
| BL layers | 18 | |
| BL thickness | 1.51 mm | **capped by the slot gap**, not δ99 (7.35 mm) |
| `--near-size` | 2.0 mm | isotropic, within 100 mm of the surface |
| wake / ground cells | 6.0 / 2.5 mm | |
| domain | −8c … +16c, +10c tall | moving ground wall at `v_ref` |

### The slot-gap cap on the boundary layer

Total BL thickness is capped at **40% of the measured minimum slot gap**, so the boundary layers extruded from *both* sides of a slot cannot collide inside it. The region between the BL cap and δ99 is covered by the isotropic slot refinement.

This cap is load-bearing, and it was being violated. `n_layers` was computed with `ceil()` and the thickness then recomputed from the geometric series — which **overshoots**: cap 1.72 mm, actual 1.83 mm. Two opposing extrusions consumed **85%** of the 4.3 mm flap1→flap2 slot, leaving 0.64 mm of inviscid core where 20% of the gap was intended. Now `floor()`: 18 layers, 1.51 mm, 70% of the slot, **1.27 mm of core**. A cap that is exceeded by construction is not a cap.

If the slot is tight enough to force fewer than 8 layers, the mesher now says so loudly — that mesh is not fit for transition-SST, and you should open the gap, raise `--yplus`, or knowingly accept wall functions.

## Which file goes into which solver

| Solver | File | Import path |
|---|---|---|
| SU2 | `.su2` | native (`MESH_FILENAME`), markers already named |
| Ansys Fluent | `.unv` | File → Import → I-deas Universal (2D solver) |
| Ansys Fluent | `.cgns` | File → Import → CGNS (alternative if UNV misbehaves) |
| OpenFOAM | `.msh` | `gmshToFoam` (file is written as msh v2.2 for this) |
| ParaView (inspect) | `.vtk` | direct open |

Boundary zones in every format: `inlet, outlet, top, ground, main, flap1,
flap2`, interior `fluid`. Run the ground as a MOVING WALL at freestream
velocity, `top` as symmetry/slip.

## Rigging conventions (the optimizer's variables)

Specified in the lift-up airfoil frame (the frame the screen ranked in);
the tool mirrors everything into car orientation. Positive = more load.

| Parameter | Default | Meaning |
|---|---|---|
| `--alpha-main` | 3.0 | main incidence [deg] |
| `--delta1` / `--delta2` | 20 / 25 | flap deflection vs PARENT chord [deg] |
| `--gap1` / `--gap2` | 0.025 | TRUE min slot gap, fraction of parent chord |
| `--overlap1` / `--overlap2` | 0.030 | flap LE upstream of parent TE, fraction of parent chord |
| `--ride-height` | 0.050 | lowest point above ground [m] |

`gap` is enforced as the measured minimum surface-to-surface distance
(bisection on the flap position), not a placement offset — the number in the
report is the number a slot gauge would read. Defaults are conventional FSAE
starting values: the CENTER of the optimization, not an answer.

## Boundary-layer sizing

First-cell height from `--yplus` (default 1.0) at `--v-ref` (default 27 m/s,
V_hi = thinnest BL, conservative), flat-plate correlation on the main chord.
y+~1 exists because the solver should run a TRANSITION model (Fluent:
Transition SST / gamma-Re_theta; SU2: SST with transition): the entire screen
logic is built on laminar/transition physics at Re 120k–510k, and fully
turbulent wall functions would throw that away.

BL total thickness is capped at 40% of the smallest measured slot gap so
opposing extrusions cannot collide inside the slot; isotropic slot
refinement covers the rest. Quality note: judge the mesh on the TRIANGLE
SICN stats in the report — BL quads are intentionally anisotropic and score
near zero on SICN by construction.

## License sanity

The team runs full Ansys, so no cell cap applies. Default settings produce
~70k cells / ~115k cells+nodes anyway (sized by physics, not license). The
script's 512k warning only matters if someone runs this on a Student seat;
`--near-size` is the first knob to coarsen in that case.

## Known limits

- 2D: no endplates, tips, or wheel wake. Treat results as section rigging +
  ground-effect trends; validate the winner in 3D.
- Center (main-only) section deferred until the 2D nose+undertray study
  defines its environment.
- Geometry prefers `kulfan_params.json` (exact ranked sections); falls back
  to `shortlist_dat/*.dat` with a warning (coordinate round-trip). Re-run
  `airfoil_screen.py` to regenerate the json if it is missing.
