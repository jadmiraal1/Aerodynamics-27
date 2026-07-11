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

## 5. Run

```
python airfoil_screen.py --quick     # smoke test first (~15 s)
python airfoil_screen.py             # full screen (~1-2 min)
```

Results land in `screen_results/` (spreadsheet, plots, XFLR5 .dat files).

## Useful flags

| Flag | Effect |
|---|---|
| `--out DIR` | write results elsewhere (e.g. `--out runs/2026-07-10`) |
| `--quick` | 150-airfoil subset, fast sanity check |
| `--no-variants` | skip CST variant generation |
| `--model-size xlarge` | more accurate NeuralFoil model (slower) |
| `--seed N` | change the variant-generation random seed |

## Everyday workflow

1. Edit the CONFIG block / `PROFILES` dict at the top of `airfoil_screen.py`
   (speeds, chords, CL targets, weights).
2. Activate the venv, run the script.
3. Inspect `screen_results/results.xlsx` and the plots.
4. Commit what you want to keep. Consider `--out` with a dated folder if you
   want to keep multiple runs side by side.

## Notes

- The `.venv` folder should NOT be committed (it's in `.gitignore`).
- No internet needed at runtime: the UIUC airfoil database ships inside the
  `aerosandbox` package.
- If `pip install` fails building numpy/scipy on an old Python, upgrade pip
  first: `python -m pip install --upgrade pip`.
