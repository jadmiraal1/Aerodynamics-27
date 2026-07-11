# FSAE Front-Wing Airfoil Screening Pipeline — Per-Element Edition

NeuralFoil-based bulk screen of the UIUC database (~2,175 airfoils) plus CST/Kulfan-perturbed variants, ranked separately per element and span station (center vs. outboard).

## Setup

```
pip install neuralfoil aerosandbox pandas openpyxl matplotlib tqdm
python airfoil_screen.py
```

Runs in ~1–2 minutes. `--quick` for a 150-airfoil smoke test, `--no-variants` to skip variant generation, `--out DIR` to change output location.

## Where the numbers come from

Speeds are from the TR26 lap sim (QSS point-mass, `vehicle_params.m` values): aero-weighted (time × v²) speed distribution on the digitized autocross/endurance tracks is 15–25 m/s with median ~19; padded for this year's expected pace to a design band of **16 / 20 / 27 m/s** (lo / design / hi).

Chords are from the packaging study (12.5 in bulkhead-to-aero-boundary, bulkhead floor 4.8 in off ground, front wheel 21.5 in behind bulkhead), using a conventional ~55/30/20 three-element split with stagger running back along the nose:

| Profile | Chord | Re band (16–27 m/s) | CL target | Role |
|---|---|---|---|---|
| Main center | 275 mm | 300k–510k | 1.10 | Center span, main plane only: lightly loaded, wake quality → undertray inlet; strict Cm |
| Main outboard | 275 mm | 300k–510k | 1.50 | Wake quality + aero balance under the flap stack; spar thickness (t/c 8–16%) |
| Flap 1 | 170 mm | 185k–315k | 1.80 | High usable CL in main element's field; Cm relaxed (outboard only) |
| Flap 2 | 110 mm | 120k–205k | 1.75 | CLmax-dominated; thin sections allowed (t/c 4.5–12%) (outboard only) |

Concept assumption: flaps run outboard only; the center span is main plane alone feeding the undertray. If `main_center` and `main_outboard` rankings converge on the same section, one mold serves the full span with incidence/chord tailoring; if they diverge strongly, that quantifies the value of a section change or dropped center plane.

Each element has its own Re list, CLmax gate, t/c gates, CL target, and scoring weights — all in the `PROFILES` dict at the top of the script.

## Pipeline stages

1. **Geometry gates** — union t/c range 4.5–16%, camber ≥ 1.5%.
2. **NeuralFoil screen** — polars at all unique Re (9 values), α −2°→18°. Per element, an airfoil must pass that element's t/c range, have confident polars at that element's Re list, and meet its CLmax gate. An airfoil can qualify for one element and not another.
3. **Variant generation** — top 6 seeds per element (deduped), 30 CST perturbations each, re-screened.
4. **Ranking** — normalized weighted score per element, computed only over that element's qualifiers.

## Metric definitions

- `*_CLmax` — worst-case CLmax across the element's Re band (low-Re end usually governs).
- `*_CL_usable` — CL at 2° below stall: usable lift with margin.
- `*_stall_gentle` — CL retained 3° past stall (1.0 = flat-top).
- `*_Cm` — pitching moment at the usable point.
- `*_LD_at_CL` — worst-case L/D at the element's target CL; zero if the section can't reach it.

## Outputs

- `results.xlsx` — `all_survivors` + `main_top` / `flap1_top` / `flap2_top` + `config`.
- `plots/<name>__<element>.png` — polars at that element's Re band.
- `shortlist_dat/*.dat` — Selig format, XFLR5-ready.

## Important caveats

- **NeuralFoil is a filter, not a verdict.** Validate the shortlist in XFLR5/XFoil at the same Re, then CFD.
- **Flap sections are screened as isolated elements.** Real flaps sit in the main element's circulation field with slot-refreshed boundary layers — they'll sustain more loading than isolated polars suggest. The profile weighting (CLmax-heavy, Cm-relaxed) approximates this, but slot/gap/overlap design needs multi-element analysis (MSES-class or CFD; XFLR5 can't do multi-element).
- **Spanwise tailoring** (center feeding the undertray vs. outboard outwash) is best handled with chord/incidence/flap-span distribution of the *same* sections, not different airfoils per span station — 