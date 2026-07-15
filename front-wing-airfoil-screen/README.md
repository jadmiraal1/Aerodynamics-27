# FSAE Front-Wing Airfoil Screening Pipeline — Per-Element Edition

NeuralFoil-based bulk screen of the UIUC database (~2,175 airfoils) plus CST/Kulfan-perturbed variants, ranked separately per element and span station (center vs. outboard). Candidates that only work in a clean wind tunnel are disqualified (see **The abuse gate**).

> **Why is anything set the way it is?** See **[`DESIGN_JUSTIFICATION.md`](DESIGN_JUSTIFICATION.md)** — every number, tagged *derived / tested / judgment / assumed*, with sensitivities and known weaknesses. This README covers *how*; that document covers *why*. Environment setup is in [`SETUP.md`](SETUP.md); the CFD mesh stage is in [`MESHING.md`](MESHING.md).

## Setup

```
pip install -r requirements.txt   # includes gmsh for the CFD mesh stage

python run_pipeline.py            # runs every stage below, in order (~4-5 min)
```

Or run the stages individually:

```
python aero_targets.py                    # 0. load budget        -> front wing CL target
python airfoil_screen.py                  # 1. screen + rank      -> results.xlsx
python compare.py                         # 2. comparison sheets  -> plots/compare_*.png
python robustness.py                      # 3. transition sweep   -> plots/robustness_*.png
python weight_sensitivity.py              # 4. weight robustness  -> plots/weight_sensitivity_*.png
python validate.py --element all          # 5. XFoil check        -> validation/   (needs xfoil)
python mesh_section.py                    # 6. CFD mesh           -> mesh_out/   (see MESHING.md)
```

Beyond the core pipeline: `shape_optimize.py` (gradient airfoil shape refinement, see it and `validate.py`'s header), `optimize.py` (2-D rigging DOE, see `OPTIMIZATION.md`), `aero_targets.py` (load budget).

**Stages 2–5 cannot run until stage 1 has written `results.xlsx`, `run_config.json` and `kulfan_params.json`.** That is why `run_pipeline.py` exists and why a VS Code *compound* launch config would not work — compounds start their configurations in **parallel**, so the downstream stages would fire before the screen had produced anything to read.

```
python run_pipeline.py --quick                   # fast chain test (rankings NOT valid; skips mesh)
python run_pipeline.py --skip-mesh               # stop after stage 4
python run_pipeline.py --only 2 3 4              # re-run a subset against existing results
```

*Linux:* if `import gmsh` fails with `libGLU.so.1: cannot open shared object file`, install the system GL libs — see [`SETUP.md`](SETUP.md). Windows/macOS wheels bundle them.

`.vscode/launch.json` has a debug configuration per stage, plus **"Pipeline: full run"** and **"Pipeline: quick smoke test"**. All figures land in `<results>/plots/`.

## The scripts

| Script | Does | Reads | Writes |
|---|---|---|---|
| `aero_targets.py` | Vehicle downforce/moment budget → what the front wing must actually produce, and what the center loading costs the outboard | vehicle assumptions at top of file | *(console)* |
| `airfoil_screen.py` | Screens the UIUC pool + CST variants, gates, scores, ranks per element | UIUC database (via AeroSandbox) | `results.xlsx`, `run_config.json`, `kulfan_params.json`, `plots/`, `shortlist_dat/` |
| `compare.py` | One side-by-side sheet per element: geometry, CL–α, drag polar, score bars | the screen's output dir | `compare_<element>.png` |
| `robustness.py` | Re-runs the shortlist under three transition environments; slope charts + retention ratios | the screen's output dir | `robustness_<element>.png`, `robustness_summary.csv` |
| `weight_sensitivity.py` | Is the airfoil choice actually driven by the scoring weights? Monte-Carlo over plausible weightings → win probability per section | the screen's output dir | `weight_sensitivity_<element>.png`, `weight_sensitivity.csv` |
| `validate.py` | Runs **XFoil** at the race condition and plots it against NeuralFoil, per section. NeuralFoil is a surrogate *for* XFoil — this is the reality check. Needs the `xfoil` binary. | the screen's output dir | `validation/<element>_validation.png`, `validation_summary.csv` |
| `mesh_section.py` | Rigs shortlist sections into the outboard 3-element stack; 2-D RANS mesh in ground effect | `kulfan_params.json` + `run_config.json` (applies the TE blunt) | `mesh_out/` |
| `common.py` | Shared config, element profiles, geometry cache, polar-reduction math | — | — |

## Sections are blunted to a manufacturable edge *before* screening

The wings are hot-wire cut from foam, and a knife-edge TE cannot be cut. Our limit is **1.5 mm** (`CONFIG['te_min_mm']`).

That is an *absolute* thickness, so it depends on chord: 1.5 mm is 0.55% of main's 275 mm chord but **1.36%** of flap2's 110 mm. Real airfoils have a near-zero TE, so **every candidate must be blunted** — with a naive "chop ≤3% of chord" gate, flap2 returned **zero** qualifiers out of 1,433.

So each candidate's `TE_thickness` is set to 1.5 mm **at that element's chord, before its polars are run**. The section that gets ranked is the section that gets cut, and the aerodynamic cost of a real edge is measured rather than assumed. Nothing is disqualified for manufacturability — everything is *made* manufacturable, and sections that depend on a feather edge pay for it in the score.

This is not cosmetic: it moved `main_outboard` from `s1210` (0.73 mm thick at 95% chord — un-cuttable) to `goe63`.

`--te-min-mm` changes the limit; `--no-mfg-gate` screens the as-drawn sharp sections.

**`kulfan_params.json` stores the as-drawn geometry** — the blunt is per-element, so one per-name file cannot carry it. Every consumer (`compare.py`, `robustness.py`, `mesh_section.py`, the polar plots, the `.dat` exports) re-applies the identical blunt from `te_min_mm` in `run_config.json`. The `.dat` files are written **per element and already blunted** (`<name>__<element>.dat`), because the same airfoil is a different section at 275 mm and at 110 mm.

> **If a geometry plot ever shows a knife edge, something has stopped applying the blunt.** It's the fastest visual check that the pipeline is self-consistent — and it caught a real bug: for a while `compare.py` was plotting `s1223_v19` with its as-drawn **0.13 mm** trailing edge, when the section that had actually been screened had a **1.50 mm** one.

## L/D is scored across a CL *band*, not at a point

The wing's operating CL is set by **incidence** — a shim at the track — and the CL target itself inherits wide uncertainty from an unfinished chain (lap sim → ClA → aero balance → undertray share → front wing load → spanwise split). Scoring L/D at a single CL let that uncertainty pick the airfoil: sweeping `main_center`'s target from 0.5 to 1.3 produced a **different winning section at every step**.

So each element has a `cl_band` (±20% of nominal), and `*_LD_band` — the **mean L/D across that band** — is what gets scored. A CL inside the band that a section cannot reach scores **zero**, penalizing it exactly in proportion to how much of the operating range it forfeits.

`*_LD_flat` (min/max of L/D across the band, 1.0 = perfectly flat) is reported as a **confidence measure on the ranking**: a flat-bucket section does not care where you trim it; a peaky one is only "best" if the target happens to be right.

`common.py` is the single source of truth. All three scripts import `PROFILES`, `CONFIG`, and the polar-reduction functions from it, so the confidence thresholds and metric definitions cannot drift apart between them.

Two handoff files make the downstream scripts honest:

- **`run_config.json`** — the config the run actually used. `compare.py` and `robustness.py` load it, so a `--model-size` override in the screen reaches them. (They used to `from airfoil_screen import MODEL_SIZE`, which binds by value at import — meaning they would silently analyze a shortlist at a different fidelity than it was ranked at.)
- **`kulfan_params.json`** — the exact geometry that was ranked. Downstream scripts analyze *these*, not the `.dat` files. The `.dat` files are a lossy coordinate export for XFLR5, and re-fitting them yields a slightly different section than the one that scored.

## Where the numbers come from

Speeds are from the TR26 lap sim (QSS point-mass, `vehicle_params.m` values): aero-weighted (time × v²) speed distribution on the digitized autocross/endurance tracks is 15–25 m/s with median ~19; padded for this year's expected pace to a design band of **16 / 20 / 27 m/s** (lo / design / hi).

Chords are from the packaging study (12.5 in bulkhead-to-aero-boundary, bulkhead floor 4.8 in off ground, front wheel 21.5 in behind bulkhead), using a conventional ~55/30/20 three-element split with stagger running back along the nose:

| Profile | Chord | Re band (16–27 m/s) | CL target | Gate (CLmax at race condition) | Role |
|---|---|---|---|---|---|
| Main center | 275 mm | 300k–510k | 1.10 | ≥ 1.10 | Center span, main plane only: lightly loaded, wake quality → undertray inlet; strict Cm |
| Main outboard | 275 mm | 300k–510k | 1.50 | ≥ 1.50 | Wake quality + aero balance under the flap stack; spar thickness (t/c 8–16%) |
| Flap 1 | 170 mm | 185k–315k | 1.80 | ≥ 1.80 | High usable CL in main element's field; Cm relaxed (outboard only) |
| Flap 2 | 110 mm | 120k–205k | 1.75 | ≥ 1.75 | CLmax-dominated; thin sections allowed (t/c 4.5–12%) (outboard only) |

Concept assumption: flaps run outboard only; the center span is main plane alone feeding the undertray. If `main_center` and `main_outboard` rankings converge on the same section, one mold serves the full span with incidence/chord tailoring; if they diverge strongly, that quantifies the value of a section change or dropped center plane.

Each element has its own Re list, CLmax gate, t/c gates, CL target, and scoring weights — all in the `PROFILES` dict in `common.py`.

## The gate: race condition, not a wind tunnel

Every candidate is analyzed a second time in the condition the car **actually races in**, at the **low-Re end** of its band (the worst corner):

```
gate case:  n_crit = 6.0,  xtr = 0.30      (open-air turbulence + light LE grit)
gate:       CLmax_gate  >=  cl_target + gate_margin        (margin = 0.00)
```

`n_crit` is Drela's e^N amplification threshold. **9 is "average wind tunnel", not free air** — outdoor autocross with wind and cars ahead churning the air sits lower, and low-Re practice lands around 5–7. `xtr = 0.30` forces transition at 30% chord: the mildly dirty leading edge a wing has after 22 km of endurance, still with some laminar run.

This gate exists because the pipeline was ranking sections its own robustness tool knew were dead. `e423` topped both `flap2` and `main_outboard`, and **five of `flap1`'s top eight were e423 variants** — yet with no laminar flow, none could reach their target CL. The failure was sitting in `robustness_summary.csv`, which nothing in the pipeline read. Now it's a gate, not a report.

**Gating on the clean case would be wrong in both directions — and "clean" isn't even the generous assumption.** At these Reynolds numbers a long laminar run means a fat laminar separation bubble, so forcing transition *earlier* often *raises* CLmax. `main_center` qualifies 702 UIUC sections at n_crit 9 but **974** at n_crit 6. A clean tunnel is a different flow, not a favorable one.

`--gate-ncrit`, `--gate-xtr`, `--gate-margin` tune the gate; `--no-gate` disqualifies nobody.

### The abuse case is reported, never gated on

A fully tripped LE (`xtr = 0.05`, rain / bugs / heavy rubber) is still computed for every candidate and reported as `*_CLmax_abuse`, `*_LD_abuse`, `*_LD_retention`. It **does not disqualify**: it's a scenario the car may never run in, and excluding on it threw away sections that race perfectly well (it cut the flap pools to 3 UIUC seeds).

Use it to break ties, and to know what you're buying. Some sections clear the race gate comfortably and still cannot reach target CL in the wet — `robustness.py` names them explicitly. That's a design decision, not a bug.

**Reality check on the flap pools:** even in a pristine wind tunnel, only ~12 UIUC sections reach flap1's target. The flap fields are small because the CL targets (1.80 at Re 186k, 1.75 at Re 120k, *isolated*) sit near the ceiling of what any section does at those Reynolds numbers — not because the gate is harsh. If an element returns fewer than `TOP_N` qualifiers the screen warns you; treat it as a signal about the target, not just the gate.

## Pipeline stages

1. **Geometry gates** — union t/c range 4.5–16%, camber ≥ 1.5%.
2. **NeuralFoil screen** — clean polars at all unique Re, α −2°→18°, plus the tripped-LE case at each element's low-Re end. Per element, an airfoil must pass that element's t/c range, have confident polars at that element's Re list, meet its CLmax gate, **and clear the abuse gate**. An airfoil can qualify for one element and not another.
3. **Variant generation** — top 6 seeds per element (deduped), 30 CST perturbations each, re-screened.
4. **Ranking** — normalized weighted score per element, computed only over that element's qualifiers.

## Metric definitions

- `*_CLmax` — **worst-case** CLmax across the element's Re band (low-Re end usually governs).
- `*_CL_usable` — CL at 2° below stall: usable lift with margin.
- `*_stall_gentle` — CL retained 3° past stall (1.0 = flat-top).
- `*_Cm` — pitching moment at the usable point.
- `*_LD_at_CL` — **SCORED.** Worst-case L/D at the element's target CL, measured **at the race condition** (n_crit 6, xtr 0.30), across the whole Re band. **NaN** (not zero) if the section can't reach the target; `*_LD_unreachable` flags this. A zero would be indistinguishable from "reached the target with catastrophic drag", and any later `mean()` would quietly lie. `score()` maps NaN → zero credit explicitly, at the point of scoring.
- `*_LD_at_CL_clean` — **reported only.** The same thing on clean polars. The gap between this and `*_LD_at_CL` is precisely how much of a section's drag advantage was laminar flow it will not have on track.
- `*_CLmax_gate` — CLmax at the race condition, low-Re end. This is what the gate tests.
- `*_CLmax_abuse`, `*_LD_abuse`, `*_LD_retention` — the fully-tripped (rain) case. **Reported, never scored.** This team does not race in the wet.

### Why L/D is scored at the race condition, not clean

The screen used to gate on the race condition but score L/D on *clean* polars — so a section could win an element on laminar-flow drag it would never have on track. It did: `main_outboard`'s top four picks (`dae31_v05`, `dae31_v14`, `mh114_v06`, `dae21_v21`) kept only **31–36%** of their L/D once the LE was dirty, while `goe63` variants ranked below them kept **60–69%**.

Scoring L/D where the car actually races removes this at the root, and needs no robustness bonus term — because there is no longer a clean-case advantage to award. After the change, `main_outboard` reorders to the `goe63` family. Absolute L/D numbers drop (≈106 → ≈87) simply because they are now quoted in the dirtier condition; that is not a regression, it is the number you were always going to get.

CLmax, `CL_usable`, `stall_gentle` and `Cm` are still taken from the clean polars, so they stay comparable with XFLR5 and published data. Only the drag-based terms moved, because drag is where laminar flow does its lying.

**Careful:** `robustness_summary.csv`'s `CLmax_clean` is measured at **one** Re (the low-Re end, or design Re with `--re design`). The xlsx's `*_CLmax` is the **worst case across the whole band**. Different quantities — don't diff them.

## Outputs

- `results.xlsx` — `all_survivors` + `<element>_top` per element + `config`.
- `run_config.json`, `kulfan_params.json` — the handoff files described above.
- `plots/<name>__<element>.png` — polars at that element's Re band.
- `shortlist_dat/*.dat` — Selig format, XFLR5-ready.
- `compare_<element>.png` — side-by-side sheets (from `compare.py`).
- `robustness_<element>.png`, `robustness_summary.csv` — transition sweep (from `robustness.py`).

## Important caveats

- **NeuralFoil is a filter, not a verdict.** Validate the shortlist in XFLR5/XFoil at the same Re, then CFD.
- **Ground effect is not modeled anywhere in this pipeline.** For a front wing at h/c ≈ 0.1–0.3 that is a real omission, and the ranking can genuinely reorder under GE. Even a shortlist that survives the abuse gate is a filter, not an answer — do not let a clean result talk you out of running the CFD.
- **Flap sections are screened as isolated elements.** Real flaps sit in the main element's circulation field with slot-refreshed boundary layers — they'll sustain more loading than isolated polars suggest. The profile weighting (CLmax-heavy, Cm-relaxed) approximates this, but slot/gap/overlap design needs multi-element analysis (MSES-class or CFD; XFLR5 can't do multi-element).
- **Post-stall behavior is the surrogate's weakest regime.** `stall_gentle` carries 20–25% of the score but rests on NeuralFoil's post-stall CL, which is smoothed. Treat it as a soft preference, not a measurement.
- **The variant pool couples the elements.** `build_variant_pool` seeds from the top-N of *every* element into one shared pool, so changing one element's profile changes which variants exist and can shift another element's ranking. Screen runs are not a controlled A/B across elements.
- **Spanwise tailoring** (center feeding the undertray vs. outboard outwash) is best handled with chord/incidence/flap-span distribution of the *same* sections, not different airfoils per span station.
