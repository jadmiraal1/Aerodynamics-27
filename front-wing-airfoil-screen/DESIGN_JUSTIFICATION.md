# Design Justification — Front Wing Airfoil Selection

**Purpose of this document.** Every number in this pipeline, where it came from, how confident we are in it, and what would change it. `README.md` explains *how* the tools work; this explains *why* they are set the way they are.

It is written to be defensible under questioning. Design judges do not ask "what does the script do" — they ask *"why 1.10?"*, *"why n_crit 6?"*, *"do you trust NeuralFoil at Re 120k?"*, and *"where is ground effect?"* Every one of those has an honest answer below, including the ones where the honest answer is "that is an assumption, and here is its sensitivity."

**The single most important thing in this document is the provenance tag on each number.** An assumption you have labelled, bounded, and sensitivity-tested is a defensible engineering decision. The same assumption presented as a derivation is a hole in your design.

---

## Provenance key

| Tag | Meaning |
|---|---|
| **DERIVED** | Computed from a model or measurement. Traceable. Changing an input changes it. |
| **TESTED** | We ran an experiment inside this repo and the number is the result. Evidence is in this doc. |
| **JUDGMENT** | An engineering choice. Not arbitrary, but not derived either. Defensible by argument, and sensitivity-bounded where it matters. |
| **ASSUMED** | A placeholder standing in for information we do not yet have. Flagged as a risk. These are the ones to attack first. |

---

## Pipeline map

| Stage | File | Output |
|---|---|---|
| Vehicle load budget | `aero_targets.py` | Front wing load target, center/outboard CL split |
| Section screen | `airfoil_screen.py` | Ranked shortlist per element (`results.xlsx`) |
| Visual comparison | `compare.py` | `compare_<element>.png` |
| Transition robustness | `robustness.py` | `robustness_summary.csv` |
| Shared config & physics | `common.py` | *(imported by all of the above)* |
| 2-D multi-element mesh | `mesh_section.py` | Solver-ready RANS mesh in ground effect (see `MESHING.md`) |

Environment setup is in `SETUP.md`.

---

## 1. Operating point

### Speeds: 16 / 20 / 27 m/s — **DERIVED**

From the TR26 lap sim (QSS point-mass). The aero-weighted (time × v²) speed distribution across the digitized autocross and endurance tracks runs 15–25 m/s with a median near 19 m/s. Padded to 16/20/27 for this year's expected pace.

The v² weighting is the right one: aerodynamic work scales with v², so the speeds where the car spends time at *high* speed matter disproportionately more than a raw time histogram would suggest.

*Weakness:* the sim is a point mass. It is adequate for a speed distribution — that is a longitudinal/grip-limited question — but see §7 for why it cannot produce an aero balance.

### Chords: 275 / 170 / 110 mm — **DERIVED**

From the packaging study: 12.5 in bulkhead-to-aero-boundary, bulkhead floor 4.8 in off ground, front wheel 21.5 in behind the bulkhead. A conventional ~55/30/20 three-element split with stagger running back along the nose.

### Reynolds bands — **DERIVED**

Re = V·c/ν with ν = 1.46 × 10⁻⁵ m²/s. Falls directly out of chord × speed:

| Element | Chord | Re band |
|---|---|---|
| main_center / main_outboard | 275 mm | 300k – 510k |
| flap1 | 170 mm | 185k – 315k |
| flap2 | 110 mm | 120k – 205k |

**Everything is evaluated worst-case across the band, not at the design point.** The low-Re end governs — it is where laminar separation bubbles are worst — and a wing that only works at 27 m/s is useless on an autocross course.

---

## 2. The concept assumption

### Center span = main plane only; flaps outboard only — **JUDGMENT**

The center span feeds the undertray, so it runs the main plane alone; the three-element stack is outboard only.

**How to test it:** if `main_center` and `main_outboard` converge on the same section, one mould serves the full span with incidence and chord tailoring. If they diverge, that quantifies the value of a section change — or of dropping the center plane entirely.

*Current result: they diverge.* `main_center` selects the `goe304` family; `main_outboard` selects `s1210`/`goe63`. That is a live design question, not a formality, and it is worth raising in design review rather than hiding.

---

## 3. CL targets — the most attackable numbers in the pipeline

### Honest history

`main_center`'s target was **1.10**, justified by a code comment reading *"lightly loaded: protect the undertray inlet."* It had **no derivation.** This mattered enormously, because we tested what happens when it changes:

**TESTED — best `main_center` section by L/D vs. the CL target you assume:**

| CL target | Winning section | Camber | L/D |
|---|---|---|---|
| 0.5 | mh43 | 0.018 | 44.9 |
| 0.7 | s7012 | 0.020 | 60.3 |
| 0.9 | e176 | 0.033 | 73.1 |
| 1.1 | fx63100 | 0.043 | 81.5 |
| 1.3 | e385 | 0.057 | 90.3 |

**A different airfoil at every step, with the winner's camber tracking the target almost linearly.** An undocumented constant was choosing the wing. This is exactly the kind of thing a judge will find, and the only defence is to have found it first.

### What the load budget actually says — `aero_targets.py`

Running the moment balance with ClA = 4.0, front axle share β = 0.40, undertray share 45%, at 20 m/s (980 N total downforce):

| Device | Front-axle load per newton produced |
|---|---|
| Front wing | **+1.47** (and −0.47 on the rear — it *lifts* the rear axle) |
| Undertray | +0.42 |
| Rear wing | **−0.19** (it *unloads* the front axle) |

**The front wing therefore needs 12–26% of total downforce, not 35–45%.** The front-axle *share* and the front-wing *share* are different quantities, because the front wing acts through a long moment arm ahead of the axle while the rear wing works against it. Sizing the wing to make 40% of total downforce would have built it roughly twice too large and pushed the car into deep understeer at speed.

At the nominal point: front wing = 19.1% = **187 N**, ClA_fw = 0.76 m², CL_fw ≈ **1.74** (3-D wing CL over S ≈ 0.44 m²).

### Why this *still* doesn't hand us `cl_target`

Two reasons, and both should be said out loud before a judge says them:

1. **The spanwise split is a free variable.** The budget fixes the wing's *area-averaged* CL. How it is distributed between center and outboard is a design choice. And it reads backwards from intuition: **the outboard's maximum achievable CL sets the center's *minimum* loading.** If the outboard tops out near CL 2.1 in ground effect, a center at CL 0.8 is unaffordable — the outboard cannot cover the deficit. That, and not the inlet wake, is what bounds the center from below.

2. **CL_fw is a 3-D wing CL; `cl_target` is a 2-D section CL.** In ground effect and inside a multi-element stack, these are not the same number. Mapping between them requires CFD — which is exactly what `mesh_section.py` exists to enable.

### The resolution: score a **band**, not a point — **JUDGMENT**

Since the target is uncertain and the wing's operating CL is set by **incidence** (a shim at the track), we do not need a section that is optimal at one CL. We need one whose **drag bucket is wide and flat across the plausible range**.

So L/D is scored as the **mean across an operating band** (`cl_band`, ±20% of nominal), with any CL inside the band that the section cannot reach scoring **zero**. This penalizes a section exactly in proportion to how much of the operating range it forfeits, and it makes the ranking **robust to the number we do not yet have.**

`*_LD_flat` (min/max of L/D across the band, 1.0 = perfectly flat) is reported as a *confidence measure on the ranking itself*: a section with a flat bucket does not care where you trim it; a peaky one is only "best" if the CL target is right.

| Element | `cl_target` (nominal) | `cl_band` (scored) | Provenance |
|---|---|---|---|
| main_center | 1.10 | 0.90 – 1.35 | JUDGMENT, bounded by `aero_targets.py` |
| main_outboard | 1.50 | 1.25 – 1.75 | JUDGMENT |
| flap1 | 1.80 | 1.55 – 2.00 | JUDGMENT |
| flap2 | 1.75 | 1.50 – 1.95 | JUDGMENT (see below) |

**`flap2` = 1.75, previously 1.90 — TESTED.** At 1.90, only **15%** of flap2's qualifiers could reach the target at all; the other 85% scored a hard zero on L/D. The survivors had a median CLmax margin of just **+0.034** above target, so L/D was being interpolated at the stall knee where CD is nearly vertical. The result was numerical noise, not discrimination: L/D scattered **26.1 → 72.9** (2.8×) across near-identical sections. At 1.75, reach rises to ~50%, the median margin becomes +0.10 (matching `flap1`'s healthy spread), and the L/D spread collapses to 1.5×.

---

## 4. The gates

### Geometry gates — **JUDGMENT**

- **t/c ranges** per element (main 8–16%, flap1 5.5–13%, flap2 4.5–12%): structural. The main plane carries the spar; flap2 at 110 mm chord can afford to be thin.
- **camber ≥ 1.5%** — **TESTED.** Symmetric sections are excluded, and we verified this rather than assuming it:

| Section | Camber | CLmax | α_stall | L/D @ CL 1.1 (race) | \|Cm\| |
|---|---|---|---|---|---|
| naca0012 | 0.000 | 1.19 | 14.2° | **41.8** | 0.026 |
| naca0015 | 0.000 | 1.25 | 16.3° | 44.5 | 0.040 |
| goe304 | 0.063 | 1.70 | 16.7° | **78.2** | 0.070 |

**A symmetric section gives the undertray a *dirtier* inlet, not a cleaner one.** To make CL 1.1 it must sit at ~14° — essentially at stall, with CLmax 1.19 against a target of 1.10 — and it pays **41.8 L/D against 78.2**. Wake thickness scales with drag, so halving L/D at the same lift roughly doubles the momentum deficit dumped into the inlet. It would also fail our own CLmax gate of 1.30.

The symmetric section wins exactly one thing: |Cm| (0.026 vs 0.070). That is real for pitch sensitivity, but it carries 0.20 weight against L/D's 0.40.

**The correct lever for a clean inlet is *loading*, not camber** — and the correct section spec is *cambered, front-loaded, low |Cm|*, which is what the `Cm_low` weight already selects for. See §6.

### Manufacturability: every section is blunted before it is screened — **TESTED**

The wings are **hot-wire cut from foam**. A knife-edge trailing edge cannot be cut: the wire melts straight through it and the foam crumbles. Our shop limit is **1.5 mm**.

This is not a minor filter. It is a hard constraint that **real airfoils systematically violate**, because a manufacturable edge is an *absolute* thickness while an airfoil is defined in *fractions of chord*:

| Element | Chord | 1.5 mm as % of chord | Median chord we'd lose to reach it |
|---|---|---|---|
| main_center / outboard | 275 mm | 0.55% | 3.0% |
| flap1 | 170 mm | 0.88% | 5.0% |
| flap2 | 110 mm | **1.36%** | **8.0%** |

**TESTED: with a naive "chop no more than 3% of chord" gate, `flap2` returned ZERO qualifiers out of 1,433.** Not one airfoil in the UIUC database can be hot-wire cut at 110 mm chord under that rule. The constraint is real and it scales inversely with chord — which is why the smallest element hurts most.

**Why we blunt rather than gate or truncate.** Three options, and only one is honest:

- *Gate on manufacturability* — rank the full section, build a modified one. This ranks geometry we cannot build. Rejected.
- *Truncate the aft* — chop at the station that gives 1.5 mm. But `s1223`'s **thin cusped aft region IS its aft loading** — it is *why* it wins. Chopping 8–21% of chord off it destroys precisely the thing being selected for. Rejected.
- **Blunt via Kulfan `TE_thickness`** — keep the chord and camber line, add a small base. NeuralFoil evaluates `TE_thickness` natively, so the **aerodynamic cost of being manufacturable is measured, not assumed.** Adopted.

So **every candidate is blunted to a 1.5 mm edge at that element's chord before its polars are run.** The section that gets ranked is the section that gets cut. Nobody is disqualified for manufacturability, because everybody is *made* manufacturable first — and sections that lean on a feather edge pay for it in the score.

**This changed the answer.** `main_outboard` had been won by `s1210`, which is **0.73 mm thick at 95% chord with a 0.08 mm trailing edge** — physically un-cuttable. Forced to carry a real edge, it loses, and `goe63` takes the element:

| | main_outboard top 3 |
|---|---|
| Sharp TE (un-buildable) | s1210_v01, s1210_v02, goe63_v00 |
| **Blunted to 1.5 mm (buildable)** | **goe63_v16, goe63_v28, goe63_v03** |

`goe63` is also the family with the best L/D retention under contamination (§5). The manufacturing constraint and the robustness argument point at the same section — which is a good sign that it is the right one.

`s1223` survives blunting at both flaps and still wins them, at CLmax ≈ 2.24 with a real 1.5 mm edge. It earns its place now.

**Reported columns:** `*_te_as_drawn_mm` (the section's own TE), `*_te_applied_mm` (what we built it to), and `*_chord_loss_pct` — how much chord we *would* have chopped had we truncated instead. **A high `chord_loss_pct` is a warning label**: that section leans hard on a thin aft region, so its blunted polars diverge most from published data, and it is the one to validate first in XFoil.

**Everything you look at is the blunted section.** `kulfan_params.json` necessarily stores the **as-drawn** geometry — the blunt is per-element, so one per-name file cannot carry it — and for a while every *consumer* of that file was plotting and exporting the raw knife edge. The gap was real and visible:

| `s1223_v19` at flap2 (110 mm chord) | TE gap | t @ 95% chord |
|---|---|---|
| As drawn (stored) | **0.13 mm** | 0.12 mm |
| **As screened / cut** | **1.50 mm** | 1.42 mm |

So `compare.py`, `robustness.py`, `mesh_section.py`, the polar plots and the `.dat` exports **all now re-apply the identical per-element blunt** from `te_min_mm` in `run_config.json`. The `.dat` files are written per element (`<name>__<element>.dat`) and are already blunted, because the same airfoil is a different section at 275 mm and at 110 mm — one file per name could not represent both.

*If a geometry plot ever shows a knife edge again, something has stopped applying the blunt. It is the fastest visual check that the pipeline is self-consistent.*

### The race-condition gate — **JUDGMENT, evidence-backed**

```
gate case:  n_crit = 6.0,  xtr = 0.30    (open-air turbulence + light LE grit)
gate:       CLmax_gate >= cl_target + 0.00      at the LOW-Re end
```

**Why not the clean case (n_crit 9)?** Two reasons.

First, n_crit 9 is Drela's *"average wind tunnel"* — **it is not free air.** Outdoor autocross with wind, and with cars ahead churning the air, sits lower; low-Re practice lands around 5–7. `xtr = 0.30` forces transition at 30% chord: the mildly dirty leading edge a wing has after 22 km of endurance, still retaining some laminar run.

Second — and this is the counterintuitive one — **TESTED: the clean case is not even the *generous* assumption.** At these Reynolds numbers a long laminar run means a fat laminar separation bubble, so forcing transition *earlier* often *raises* CLmax:

| Transition assumption | main_center qualifiers (of 1,433 UIUC) |
|---|---|
| clean, n_crit 9 | 702 |
| **track, n_crit 6** | **974** |
| windy, n_crit 4 | 993 |

A clean tunnel is a *different* flow, not a favourable one.

**Why the gate exists at all — TESTED.** The pipeline was ranking sections its own robustness tool knew were dead. `e423` topped both `flap2` and `main_outboard`, and **five of `flap1`'s top eight were e423 variants** — yet with no laminar flow, *none of them could reach their target CL at all*. The failure was sitting in `robustness_summary.csv`, which nothing in the pipeline read. It is now a gate, not a report.

### The abuse case is reported, never gated — **JUDGMENT**

A fully tripped LE (`xtr = 0.05` — rain, bugs, heavy rubber) is computed for every candidate and reported as `*_CLmax_abuse`, `*_LD_abuse`, `*_LD_retention`. **It does not disqualify.**

Rationale: **this team does not race in the wet.** Gating on it was tested and cut the flap pools to **3 UIUC seeds** — it was throwing away sections that race perfectly well, to protect against a scenario that cannot occur. Use it to break ties, and to know what you are buying.

---

## 5. L/D is scored where the car races — **TESTED**

The screen originally *gated* on the race condition but *scored* L/D on **clean** polars. That let a section win an element on laminar-flow drag it would never have on track — and it did:

| | Section | L/D (as scored) | L/D retention when tripped |
|---|---|---|---|
| **Old** (clean-scored) | dae31_v05 | 106.4 | 0.36 |
| | dae31_v14 | 109.5 | 0.32 |
| | dae21_v21 | 110.6 | **0.31** |
| **New** (race-scored) | **goe63 / s1210 family** | ~87 | **0.64 – 0.69** |

`main_outboard` reorders completely. Absolute L/D drops (≈106 → ≈87) purely because it is now quoted in the dirty condition — **that is not a regression, it is the number we were always going to get.**

`*_LD_at_CL_clean` is retained as a reported column: **the gap between it and the scored value is precisely how much of a section's drag advantage was laminar flow it will not have.**

Note this needed **no robustness bonus term**: once you score where you race, there is no clean-case advantage left to award.

CLmax, `CL_usable`, `stall_gentle` and `Cm` are still taken from clean polars, deliberately, so they remain comparable with XFLR5 and published data. Only the drag terms moved — because drag is where laminar flow does its lying.

---

## 6. Scoring weights — **JUDGMENT**

Weights live in `PROFILES` (`common.py`). They encode each element's *job*, and they are the most explicitly subjective part of the pipeline. State them as choices, not facts.

| Element | Dominant weights | Why |
|---|---|---|
| **main_center** | LD_band 0.40, stall_gentle 0.25, **Cm_low 0.20** | Feeds the undertray. Wake quality first. Lightly loaded, so CLmax barely matters (CL_usable 0.05). |
| **main_outboard** | LD_band 0.35, stall_gentle 0.25, Cm_low 0.15 | Sets the wake and the car's aero balance; must be thick enough for a spar. |
| **flap1** | CLmax 0.35, CL_usable 0.25 | Loaded by the main element's circulation. Cm penalty relaxed — flap Cm reacts through the stack. |
| **flap2** | CLmax 0.40, CL_usable 0.25 | Smallest chord, lowest Re, most aggressive loading. Thin sections fine. |

### Why `Cm_low` on the center element is the *right* answer to "does a cambered wing dump flow into the nose?"

Flow turning is set by **circulation — i.e. CL — not by geometric angle of attack.** A section making CL 1.1 turns the flow by the same amount whether it got there via camber at 2° or symmetry at 14°. You cannot generate downforce without deflecting air upward; that deflection *is* the downforce. No geometry escapes it.

What you *can* control is **where the load sits along the chord.** A front-loaded section builds circulation near the leading edge and leaves the trailing edge relatively axial. An aft-loaded one hooks the flow hard at the back — right where it is pointed at the nose. **The signature of that is |Cm|.** So `Cm_low` at 0.20 is not a pitch-sensitivity term alone; it is the wake-direction term. If flow into the nose is the concern, that weight is the knob — not symmetry.

### NeuralFoil model size: `xxxlarge` — **TESTED**

NeuralFoil ships **8 networks**: `xxsmall`, `xsmall`, `small`, `medium`, `large`, `xlarge`, `xxlarge`, `xxxlarge`. The library's own default is `xlarge`.

This pipeline was running **`large` — below the library default — with no justification.** It was simply the value in the original script.

**The cost argument for staying small does not exist.** Full 1,433-airfoil screen:

| Model | Full run |
|---|---|
| large *(was)* | 0.6 min |
| xlarge *(library default)* | 0.7 min |
| xxlarge | 1.2 min |
| **xxxlarge** *(now)* | **3.3 min** |

2.7 minutes, once per design cycle, to use the most accurate model available.

**And it changes the answer** — in the same pattern as the weights:

| Element | Spearman ρ vs `large` | Winning family across model sizes |
|---|---|---|
| main_center | **0.61 – 0.78** | goe319 / cootie / goe304 — **flips** |
| main_outboard | 0.75 – 0.87 | s1210 — stable |
| flap1 | 0.98 | s1223 — stable |
| flap2 | 0.96 – 0.99 | s1223 — stable |

**The important part is *why* it changes.** CLmax differences between models are tiny (0.006–0.04) — the physics barely moves. What moves is the **ranking**, because the top candidates sit within ~0.01 of each other in score, so model noise exceeds the gaps between them.

**This is the same conclusion the weight study reached, from a completely different direction: the shortlist is a cloud of near-identical siblings, and "rank 1" among them is substantially noise.** Trust the *family*, not the variant. The flaps are insensitive to both the weights and the model; `main_center` is sensitive to both, and is the element to treat with the least confidence.

### Are the weights choosing the wing? — **TESTED** (`weight_sensitivity.py`)

The weights are the most subjective thing in the pipeline. That is only a *problem* if the **answer** depends on them. So we measured it: draw thousands of weight vectors from a Dirichlet distribution centred on the nominal weights, re-score, and record who wins.

**The number that matters is the FAMILY win probability, not the variant.** The rivals are almost always CST siblings of the winner (`goe63_v16` vs `goe63_v28` vs `goe63_v03`). You do not cut a variant — you cut a **family**: one mould, one section. If the weights only shuffle the order *within* a family, they are not choosing the wing; the aerodynamic decision was already made, and they are picking between near-identical siblings.

**Nominal spread (concentration 40 — weightings a reasonable engineer might have chosen):**

| Element | Winner | Variant win % | **Family win %** | Margin over 2nd | Spearman ρ | Verdict |
|---|---|---|---|---|---|---|
| main_center | goe304 | 69% | **69%** (5 families) | 0.0091 | 0.937 | MODERATE |
| main_outboard | goe63_v16 | 73% | **78%** (3 families) | 0.0180 | 0.872 | ROBUST |
| flap1 | s1223_v19 | 56% | **100%** (1 family) | 0.0059 | 0.974 | ROBUST |
| flap2 | s1223_v19 | 100% | **100%** (1 family) | 0.0166 | 0.990 | ROBUST |

**`flap1` is the clearest demonstration of why the family view is the right one.** At variant level it looks shaky — the winner takes only 56% of draws. But **every single draw is won by an `s1223` variant.** The weights are not choosing the airfoil at all; they are choosing which perturbation of `s1223` sits on top. For a mould, that is a 100% robust answer wearing a 56% disguise.

**Stress test (concentration 8 — wild disagreement about what matters):**

| Element | Family win % | Verdict |
|---|---|---|
| main_center | 34% | FRAGILE |
| main_outboard | 48% | FRAGILE |
| flap1 | 98% | ROBUST |
| flap2 | 100% | ROBUST |

**Honest reading.** The flaps are genuinely weight-insensitive — `s1223` wins under essentially any weighting, because it is simply the highest-CLmax section that survives the gates. **The mains are not.** Under wide disagreement about the weights, `main_center` and `main_outboard` both become coin flips among 7–8 families.

That is the correct and expected result, and it should be *stated*, not hidden: the flaps are a CLmax problem with one dominant answer, while the mains are a genuine multi-objective trade (wake quality vs. stall vs. Cm) where the weighting reflects a real engineering position. **`main_center` is the element whose answer most depends on judgment** — its margin over 2nd place is 0.0091, which is small.

The defensible claim is therefore: *"`s1223` is the flap answer under any weighting. For the mains, `goe63`/`goe304` win under our stated priorities, and we can show exactly how much that conclusion leans on them."* That is a much stronger position than pretending the weights are objective.

**How to use the figures.** Each `weight_sensitivity_<element>.png` has four panels: win probability (who wins across plausible weightings), one-at-a-time sweeps (which single weight, if any, the result hinges on — the winner typically holds rank 1 across a broad plateau and only loses when one metric exceeds ~0.5–0.6), the nominal ranking with its margin, and rank-order stability. The pure-criterion winners along the bottom are the corners of the weight simplex: if a section wins at several corners, no interior weighting can dislodge it.

### Scoring mechanics

Each metric is min-max normalized over the **5th–95th percentile** of that element's qualifiers, then weighted and summed. Percentile clipping (rather than min-max over the full range) prevents a single outlier from compressing everyone else into a narrow band. The cost is saturation at the top: several sections can tie at 1.0 on a metric. This is a known limitation, and it is why `*_LD_flat` matters — it distinguishes candidates the score cannot.

---

## 7. Known weaknesses — read this before design review

These are real, and they are better raised by you than found by a judge.

1. **Ground effect is not modeled in the screen.** A front wing runs at h/c ≈ 0.1–0.3, where CLmax, the drag polar, and the CoP all differ substantially from free air. **Every ranking in `results.xlsx` is a free-air ranking.** This is the single biggest limitation, and it is the reason `mesh_section.py` exists. A shortlist that survives the screen is a *filter*, not a verdict — do not let a clean screen result talk you out of the CFD.

2. **Flaps are screened as isolated sections.** Real flaps sit in the main element's circulation field with slot-refreshed boundary layers, and will sustain more loading than isolated polars suggest. The profile weighting (CLmax-heavy, Cm-relaxed) *approximates* this; it does not model it. Slot/gap/overlap design is a multi-element problem — hence `mesh_section.py`.

3. **NeuralFoil is a surrogate for XFoil, not for reality.** It inherits XFoil's own limits (an integral-BL method with an e^N transition model), and at flap2's Re ≈ 120k those limits are real: laminar separation bubbles dominate and are precisely what integral methods handle worst. We mitigate with a confidence gate (`conf_min = 0.85`, `min_valid_frac = 0.6`) and by rejecting any polar whose CLmax lands on the sweep boundary — but **NeuralFoil is a filter, not a verdict.** Validate in XFoil/XFLR5, then CFD.

4. **Post-stall behavior is the surrogate's weakest regime.** `stall_gentle` carries 20–25% of the score but rests on NeuralFoil's post-stall CL, which is smoothed. Treat it as a soft preference, not a measurement.

5. **The lap sim is a point mass, so it cannot produce an aero balance.** A point mass has no axles. β = 0.40 is therefore a **TEAM TARGET, not a derivation.** Fixing this needs a two-track (4-corner) quasi-static model with **load-sensitive tires** — without load sensitivity, downforce appears free and any CLA optimum is fiction. This is the highest-leverage item on the vehicle-dynamics side.

6. **The undertray share (45%) is ASSUMED, and it is on the critical path.** Sweeping it from 35% to 55% swings the front wing target from 223 N to 151 N — about as much as the entire 0.35–0.45 balance range does. It is also *circular*: the undertray's share partly depends on the front wing's wake. Bootstrap, then iterate after CFD.

7. **Vehicle parameters are last year's car.** Mitigation: treat them as a nominal point and sweep the plausible range (`aero_targets.py --cla / --beta / --undertray`). If the target is stable across it, last year's numbers were good enough — and *that is a result you can defend*. If it moves, you have found what to pin down first.

8. **The variant pool couples the elements.** `build_variant_pool` seeds from the top-N of *every* element into one shared pool, so changing one element's profile changes which variants exist and can shift another element's ranking. Screen runs are **not** a controlled A/B across elements.

9. **Shortlist diversity is limited.** Both flap shortlists are dominated by `s1223` perturbations. Even in a pristine tunnel only ~12 UIUC sections reach flap1's target — **the flap fields are small because the CL targets sit near the ceiling of what any section does at those Reynolds numbers**, not because the gates are harsh. Worth knowing before committing CFD hours to fifteen near-identical candidates.

---

## 8. Reproducibility

- **Versions are pinned** (`neuralfoil==0.3.2`, `aerosandbox==4.2.10`) precisely so results are reproducible across machines. This was verified: an independent run on different hardware reproduced the shortlist **byte-for-byte** on seed 42.
- **`--seed 42`** controls CST variant generation.
- **`run_config.json`** records the config each run actually used; `compare.py` and `robustness.py` load it back rather than importing constants, so they cannot analyze a shortlist at a fidelity it was never screened at.
- **`kulfan_params.json`** stores the *exact* geometry that was ranked. Downstream tools (including `mesh_section.py`) consume it. The `.dat` files are a lossy XFLR5 export — re-fitting them yields a slightly different section than the one that scored.

---

## 9. Open items

1. **Two-track QSS lap sim with load-sensitive tires.** Gates every aero target, not just this one.
2. **Pin the undertray share.** Currently the largest single unknown in the load budget.
3. **CFD in ground effect, multi-element** (`mesh_section.py` → solver). This is what converts the shortlist into a decision, and what maps 3-D wing CL onto 2-D section CL.
4. **Resolve center vs. outboard.** They diverge. Decide whether that justifies two moulds, or whether the center plane should exist at all.
5. **Consider widening the flap seed pool** before spending CFD hours — see weakness 9.
