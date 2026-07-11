# FSAE Front-Wing Airfoil Screening Pipeline

NeuralFoil-based bulk screen of the UIUC database (~2,175 airfoils) plus CST/Kulfan-perturbed variants, ranked against two spanwise role profiles for a multi-element front wing.

## Setup

```
pip install neuralfoil aerosandbox pandas openpyxl matplotlib tqdm
python airfoil_screen.py
```

Full run takes ~1 minute. `--quick` for a 150-airfoil smoke test, `--no-variants` to skip variant generation, `--out DIR` to change output location.

## The two role profiles

**CENTER (undertray feed):** the center span must deliver clean, low-loss flow to the undertray, and over-loading the front wing hurts your aero balance anyway. Scored on L/D at CL=1.4 (thin wake), stall gentleness (robust to ride-height/pitch changes), low |Cm| (less pitch sensitivity), usable CL, and max L/D.

**OUTBOARD (multi-element + outwash):** drives downforce and redirects flow around the front tires; sections here see flap-induced loading and tire-wake dirt. Scored on worst-case CLmax across the Re sweep, usable CL (2° stall margin), stall gentleness, L/D at CL=1.7, and |Cm|.

All weights, CL targets, Re list, and gates are in the CONFIG block at the top of `airfoil_screen.py`. Edit Re to match your actual chord lengths and speed range (Re ≈ 68,500 × V[m/s] × c[m] at sea level).

## Pipeline stages

1. **Geometry gates** — t/c between 5.5% and 16%, camber ≥ 1.5%. Kills ~1/3 of the database.
2. **NeuralFoil screen** — polars at Re 200k/350k/500k, α −2°→18°. Points with analysis_confidence < 0.85 discarded; airfoils with CLmax on the sweep boundary or worst-case CLmax < 1.4 gated out. This eliminates the large majority.
3. **Variant generation** — top 15 seeds (combined score) fitted to Kulfan/CST parameters, 30 Gaussian perturbations each, re-screened. Variants often beat their seeds by trading a little CLmax for gentler stall / lower Cm.
4. **Ranking + export** — normalized weighted scores per role.

## Outputs

- `results.xlsx` — `all_survivors` (everything that passed gates, all metrics), `center_top` / `outboard_top` (top 15 per role), `config` (exact settings used, for traceability).
- `plots/*.png` — CL–α, drag polar, L/D for each shortlisted airfoil at all three Re.
- `shortlist_dat/*.dat` — Selig-format coordinates, importable directly into XFLR5.

## Metric definitions

- `CLmax_worst` — minimum CLmax across the Re sweep (worst case, not average).
- `CL_usable` — CL at 2° below stall: what you can actually run with margin.
- `stall_gentle` — CL retained 3° past stall (1.0 = flat top, no drop).
- `Cm_use` — pitching moment at the usable point.
- `LD_at_CL_center/outboard` — worst-case L/D at the role's target CL; zero if the airfoil can't reach that CL.

## Important caveats

- **NeuralFoil is a filter, not a verdict.** It's a surrogate trained on XFoil, so it inherits XFoil's limits (2D, free transition, no ground effect) plus its own approximation error. Validate the shortlist in XFLR5/XFoil, then CFD.
- **Single-element polars ≠ multi-element behavior.** Slot geometry, flap overlap/gap, and ground effect dominate a real front wing. Use this screen to pick candidate *sections*; the element decomposition and slot design is a separate optimization (XFLR5 can't do multi-element — you'll want MSES-style analysis or CFD).
- **S1223-class sections rank mid-pack here by design** — your criteria penalize the high Cm and sharp stall that come with maximum-camber sections. If you want a pure max-downforce ranking, sort `all_survivors` by `CLmax_worst`.
- Variant `.dat` files are analytically smooth (Kulfan) but haven't been checked for TE closure quirks — inspect before manufacturing.

## Next steps after screening

1. Import `shortlist_dat` into XFLR5, run Type-1 polars at your Re list, confirm NeuralFoil agreement (expect ±5% on CLmax at these Re).
2. Pick main-element vs. flap sections (flaps can tolerate higher camber/Cm).
3. Multi-element + ground-effect CFD on 2–3 configurations.
