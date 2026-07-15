# Environment Setup (Windows)

One-time setup to run the airfoil screen locally.

## 1. Check Python

You need Python 3.9–3.12. In PowerShell:

```
py --list
```

If nothing 3.9–3.12 shows up, install from https://www.python.org/downloads/
(check "Add python.exe to PATH" during install).

## 2. Create a virtual environment (once)

From this folder (`front-wing-airfoil-screen`):

```
py -3.12 -m venv .venv
```

(Use whatever 3.9–3.12 version you have, e.g. `py -3.11`.)

## 3. Activate it (every new terminal)

```
.venv\Scripts\Activate.ps1
```

If PowerShell refuses with an execution-policy error, run this once and retry:

```
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
```

You'll know it worked when `(.venv)` appears at the start of your prompt.

## 4. Install dependencies (once, ~1 min)

```
pip install -r requirements.txt
```

This installs everything, **including `gmsh`** for the CFD meshing stage. There is no separate install step for it.

### If you are on Linux (not needed on Windows/macOS)

The `gmsh` wheel links against system OpenGL libraries that pip cannot install. If `import gmsh` fails with `libGLU.so.1: cannot open shared object file`:

```
sudo apt-get install libglu1-mesa libxcursor1 libxinerama1 libxft2
```

Windows and macOS wheels bundle these — nothing extra to do.

## 5. Run

```
python run_pipeline.py --quick     # smoke test first (~20 s, no mesh)
python run_pipeline.py             # everything (~4-5 min)
```

`run_pipeline.py` runs all six stages **in order**. It has to: stages 2–5 read `results.xlsx`, `run_config.json` and `kulfan_params.json`, and only stage 1 produces them.

| Stage | Script | Produces |
|---|---|---|
| 0 | `aero_targets.py` | vehicle load budget → what the front wing must make |
| 1 | `airfoil_screen.py` | `results.xlsx` + the two handoff JSONs |
| 2 | `compare.py` | `plots/compare_*.png` |
| 3 | `robustness.py` | `plots/robustness_*.png`, `robustness_summary.csv` |
| 4 | `weight_sensitivity.py` | `plots/weight_sensitivity_*.png` |
| 5 | `mesh_section.py` | `mesh_out/` — 2-D RANS mesh in ground effect |

Everything lands in `screen_results/`. Figures are all in `screen_results/plots/`.

**The `--quick` run's rankings are meaningless** — it screens 150 of 2,175 airfoils. It exists to prove the chain runs, and it writes to `screen_results_smoke/` so it can never overwrite real results. It skips the mesh stage for the same reason: meshing a meaningless rank-1 section wastes a minute.

## Useful flags

| Flag | Effect |
|---|---|
| `--quick` | 150-airfoil subset, fast chain test. Implies `--skip-mesh`. |
| `--skip-mesh` | stop after stage 4 |
| `--only 2 3 4` | re-run just those stages against an existing results dir |
| `--results DIR` | write/read somewhere else (e.g. `--results runs/2026-07-13`) |
| `--no-variants` | skip CST variant generation |
| `--seed N` | change the variant-generation random seed |

Individual stages take their own flags — `python airfoil_screen.py --help`. Notably `--te-min-mm` (hot-wire limit), `--gate-ncrit` / `--gate-xtr` (racing condition), and `--model-size`.

## Everyday workflow

1. Edit **`common.py`** — `PROFILES` (chords, CL targets/bands, gates, weights) and `CONFIG` (model size, transition assumptions, TE limit). *This is the only file with knobs in it*; the other scripts import from it, so they cannot drift apart.
2. Activate the venv, `python run_pipeline.py`.
3. Read `screen_results/results.xlsx` and `screen_results/plots/`.
4. Before you believe a ranking, read **`DESIGN_JUSTIFICATION.md`** — it says which numbers are derived, which are assumed, and how much the answer leans on them.

## VS Code

`.vscode/launch.json` has a debug configuration per stage, plus **"Pipeline: full run"** and **"Pipeline: quick smoke test"**.

A VS Code *compound* will not work for the full chain — compounds launch their configurations **in parallel**, so stages 2–5 would fire before the screen had written anything. That is what `run_pipeline.py` is for.

## Notes

- The `.venv` folder should NOT be committed (it's in `.gitignore`).
- Generated output (`plots/`, `shortlist_dat/`, `mesh_out/`) is also gitignored. It is derived data — rebuild it, don't commit it. The screen **purges those directories on every run**, so stale artifacts from an older run cannot survive and be mistaken for current ones.
- No internet needed at runtime: the UIUC airfoil database ships inside the `aerosandbox` package.
- If `pip install` fails building numpy/scipy on an old Python, upgrade pip first: `python -m pip install --upgrade pip`.
