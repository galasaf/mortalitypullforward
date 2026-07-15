# Mortality Pullforward Model

## What This Model Does

This model quantifies the **mortality pullforward effect** from COVID-19: because COVID disproportionately killed people who were already in poor health (and therefore had shorter remaining life expectancy), the surviving population post-COVID is healthier on average — leading to lower-than-baseline mortality and longer life expectancy for survivors.

---

## Quick Start

### Interactive interface (no terminal needed)

Double-click **`interface.html`** to open the Mortality Pullforward Explorer in a
browser. It contains a JavaScript port of the exact same model (verified to match
the Python results to ~1e-14, including improvement and the equivalent
multiplier) with:

- A **Model mode** switch at the top of the sidebar: "Assume pullforward"
  (the original mode, everything below) or "**Calibrate from excess**" — the
  reverse direction, where the user assumes each 5-year age band's cumulative
  COVID excess deaths (% of one normal year's deaths, default 50%, all
  attributed to 2020) and the model *solves* the pullforward: the share of
  2021's deaths pulled into 2020, grading linearly to zero over a harvest
  grade-out horizon. Includes a valuation-year input (default 2025; LE gain
  and multiples are for people alive *then*), a 2010–2035 mortality-trajectory
  chart (baseline vs COVID line, improvement projected backward and forward
  from the 2019 table), an optional toggle grading the 2020 excess spike over
  N years (re-times the spike only; totals and the solve are unchanged), an
  infeasibility warning when the excess exceeds all deaths within the horizon,
  and its own tiles/charts/summary table. Its Python twin is `excess.py`
  (verified to ~1e-14 against the JS). See "Calibration mode" below.
- All 9 named scenario presets, plus free-form sliders for peak fraction,
  grade-out horizon (0–100; 0 = no pullforward), shape (linear/step/exponential),
  and decay rate. Step = full peak effect for every year inside the horizon
  (a box); if that removes 100% of a cohort's remaining deaths the UI shows a
  "no survivors" warning instead of numbers
- **Custom segmentation**: user-defined age cutoffs and an optional sex split
  define input segments (e.g. "<65 / 65+" or "8 segments: 4 age bands × M/F");
  every segment carries its own full set of pullforward inputs plus an on/off
  switch, edited via the "Editing segment" selector in the sidebar. Cohorts are
  matched to segments by their age at the end of COVID. Presets map onto this
  (age-varying presets become 4 segments at cutoffs 65/75/85, labeled
  45–64/65–74/75–84/85+)
- A left input pane that scrolls independently of the page, with collapsible
  sections, and an **Output view** bar above the results (focus cohort for
  tiles/charts; sex checkboxes and an editable age list for the summary table)
- Live-updating life expectancy tiles, charts of f(t), the mortality multiple,
  and d(t) vs g(t), and a summary table across all ages/sexes
- A **Mortality improvement** sidebar group: on/off toggle, flat %/yr rate
  (default 1%), CSV import of 1D/2D scales (same schema as the Python tool),
  and client-side template downloads
- An **Equivalent mortality multiplier** tile and summary-table column
- A generated `sensitivity.py` / `scenarios.py` command matching the current
  slider settings (including `--improvement-*` flags), so any interactive
  result can be reproduced in Python
- A built-in "How the model works" explainer of the four-step math

The Python code below remains the source of truth for saved CSV output and sweeps.

**Hosting it online:** `interface.html` is fully self-contained (no server, no
external requests), so any static host works. Two deployments exist:

1. **GitHub Pages (live)**: the repo is `github.com/galasaf/mortality-pullforward`
   and the tool is served at **https://galasaf.github.io/mortality-pullforward/**
   (the root `index.html` is just a redirect to `interface.html`). Updating it =
   commit + `git push` to `main`; Pages rebuilds automatically in ~1 minute.
2. **iWebFusion** shared hosting: `deploy/index.html` is the upload-ready copy
   (cPanel → File Manager → `public_html`); see "How to run this tool.txt"
   Part 5 for the steps. **After changing `interface.html`, refresh the copy**
   (`cp interface.html deploy/index.html`) and re-upload it.

Nothing entered into the page leaves the visitor's browser; imported CSVs are
parsed locally.

### Testing sensitivities (easiest — no file editing)

Use `sensitivity.py` to try any assumption or sweep a parameter from the
command line. **This is the recommended way to test "what if" questions.**

```bash
# See all options and worked examples
python sensitivity.py

# Single run with custom assumptions
python sensitivity.py --peak 0.6 --grade-out 7 --age 65
python sensitivity.py --shape exponential --decay 0.3

# Tweak an existing preset
python sensitivity.py --base moderate_base --peak 0.8

# SWEEP one knob across a range and see how life expectancy responds
python sensitivity.py --sweep peak      --range 0.3 0.8 0.1 --age 65
python sensitivity.py --sweep grade-out --values 2 5 7 10 15
python sensitivity.py --sweep decay     --range 0.1 0.5 0.1 --shape exponential

# MORTALITY IMPROVEMENT (default: flat 1%/yr at every age and year)
python sensitivity.py --improvement-rate 0.015           # flat 1.5%/yr instead
python sensitivity.py --no-improvement                   # static life table (old behavior)
python sensitivity.py --improvement-template 1d          # write a 1D CSV template & exit
python sensitivity.py --improvement-template 2d          # write a 2D CSV template & exit
python sensitivity.py --improvement-table my_scale.csv   # import a filled-in 1D or 2D CSV
python sensitivity.py --sweep improvement --range 0.0 0.02 0.005 --age 65
```

Knobs: `--peak` (peak_fraction), `--grade-out` (horizon years), `--decay`
(exponential rate), `--shape` (linear|step|exponential), `--improvement-rate`
(flat annual improvement). Filters: `--age`, `--sex`.
Sweeping `peak`/`grade-out`/`decay`/`improvement` prints a table of LE change
vs the swept value.

### Calibrating from COVID excess (the reverse direction)

Use `excess.py` when the starting point is *how many extra deaths COVID
caused* rather than an assumed pullforward curve:

```bash
python excess.py                       # print full help + worked examples
python excess.py --run                 # defaults: 50% excess, 7-yr grade-out, valued 2025
python excess.py --excess-all 60 --grade-out 10 --valuation-year 2027
python excess.py --excess "30,30,...,45"   # 21 values, one per 5-yr band (0-4 ... 100+)
python excess.py --grade-excess 3      # spread the 2020 spike over 3 years (timing only)
python excess.py --age 65 --sex male --save   # CSVs to output/excess_calibration/
```

It prints the solved "pulled from 2021" share per cohort (much smaller than
the excess % — see Calibration mode below), LE gain at the valuation year,
the equivalent flat multiplier, and a 2010–2035 mortality trajectory table.
Improvement flags (`--improvement-rate/-table`, `--no-improvement`) work as
in `sensitivity.py`. Ages passed with `--age` are ages **in 2020**.

### Running named scenarios

Open a terminal in this folder and run:

```bash
# Print results for the default scenario
python main.py

# Print AND save CSVs to output/default/
python main.py --save

# Run a specific named scenario
python scenarios.py moderate_base

# Run a scenario and save CSVs to output/moderate_base/
python scenarios.py moderate_base --save

# Compare all scenarios side by side (LE change table)
python scenarios.py compare_all

# Compare all scenarios, filtered to age 65 only
python scenarios.py compare_all 65

# List all available scenarios with descriptions
python scenarios.py
```

---

## Output

By default results print to the terminal only. With `--save`, CSVs are written to:

```
output/
└── <scenario_name>/
    ├── summary.csv                   ← one row per age/sex: LE pre, post, change,
    │                                    equiv_flat_qx_multiplier, survivor %
    ├── multiples_male_age65.csv      ← year-by-year: qx, pullforward f(t), d(t), g(t), multiple
    ├── multiples_female_age65.csv
    └── ...                           ← one file per age × sex combination
```

`equiv_flat_qx_multiplier` is the **equivalent flat mortality multiplier**: the
single factor `m` such that multiplying qx at *every* age of the baseline table
by `m` reproduces the survivors' post-COVID LE. E.g. `0.92` means the
survivors' LE gain is equivalent to a uniform 8% cut in mortality across the
whole table. It is solved by bisection (exact to ~1e-6 years) against the same
baseline used for `le_pre_covid` (i.e. including mortality improvement when
enabled). It also appears in every printed summary as "Equiv. Mort. Mult. (%)".

---

## How to Modify Assumptions

| What you want to change | Where |
|---|---|
| Pullforward shape, grade-out, peak fraction | `config.py` → edit any preset function's `PullforwardConfig` |
| Add a new scenario | Add a `def my_scenario_config()` in `config.py`, register it in `SCENARIOS` dict in `scenarios.py` |
| Ages / sexes to analyze | `config.py` → `ModelConfig.analysis_ages` / `.analysis_sexes` |
| Use a real SSA mortality table | `config.py` → `ModelConfig.ssa_table_path = "path/to/file.csv"` (columns: `age`, `male_qx`, `female_qx`) |
| Mortality improvement (flat rate, on/off, CSV table) | `config.py` → `ModelConfig.improvement_rate` / `.improvement_enabled` / `.improvement_table_path`, or the `--improvement-*` / `--no-improvement` flags on `sensitivity.py` |
| Pullforward math itself | `mortality_model/pullforward.py` → `compute_pullforward_fraction()` |
| Life expectancy formula | `mortality_model/analysis.py` → `compute_life_expectancy()` |
| Excess-calibration math (solved pullforward, trajectory) | `mortality_model/excess.py` → `run_excess_cohort()` / `mortality_trajectory()` |

---

## Available Scenarios

| Scenario name | Peak fraction | Grade-out | Age-varying | Description |
|---|---|---|---|---|
| `default` | 100% | 10 yrs | Yes | Theoretical ceiling; peak=100% is aggressive |
| `short_harvest` | 60% | 3 yrs | No | Conservative: flu-harvesting view |
| `moderate_base` | 65% | 7 yrs | Yes | Most defensible for COVID |
| `elderly_concentrated` | 70% | 3–15 yrs | Strong | 85+ steep, <65 flat |
| `long_harvest` | 50% | 15 yrs | Yes | Aggressive upper bound |
| `exponential_moderate` | 70% | fat tail (decay=0.4) | No | No hard cutoff |
| `exponential_long_tail` | 80% | fat tail (decay=0.2) | No | Slow decay, broad effect |
| `sickest_only` | 90% | 2 yrs | No | Only the imminently dying; nearly all of them |
| `uniform_benchmark` | 50% | 10 yrs | No | Clean benchmark, no age-variation |

---

## Key Parameters (`config.py`)

### `PullforwardConfig`

| Parameter | Default | Description |
|---|---|---|
| `peak_fraction` | `1.0` | Scales the entire f(t) curve. `peak_fraction=0.6` means 60% of year-1 deaths were pulled into COVID. Realistic range: 0.3–0.8. |
| `shape` | `'linear'` | `'linear'` = ramps from peak to zero at grade-out; `'step'` = full peak for every year inside the horizon, zero after; `'exponential'` = fat tail, never reaches zero |
| `default_grade_out_years` | `10` | For linear/step: the horizon in years. `0` disables the pullforward entirely (f≡0) |
| `exponential_decay_rate` | `0.3` | For exponential: higher = steeper decay |
| `age_varying` | `True` | If True, older cohorts use a shorter grade-out |
| `age_grade_out` | see config | Dict mapping age bands to grade-out years |

### `ModelConfig`

| Parameter | Default | Description |
|---|---|---|
| `covid_end_year` | `2022` | Reference point for the projection |
| `ssa_table_path` | `None` | Custom SSA CSV path; None uses built-in SSA 2019 approximation |
| `max_age` | `119` | Oldest age in life table |
| `analysis_ages` | `[40,50,55,60,65,70,75,80]` | Cohort starting ages to analyze |
| `analysis_sexes` | `['male','female']` | Sexes to analyze |
| `improvement_enabled` | `True` | Apply future mortality improvement to the baseline table |
| `improvement_rate` | `0.01` | Flat annual improvement rate (1%/yr at every age/year) used when no CSV is supplied |
| `improvement_table_path` | `None` | Path to a 1D or 2D improvement CSV (schema auto-detected; see below) |

---

## Mortality Improvement

The model projects future mortality improvement on top of the base life table.
**Default: flat 1%/yr improvement at every age and year.** Implementation:
`mortality_model/improvement.py`.

**Convention.** Improvement compounds from `covid_end_year` (2022) forward
along each cohort's age/calendar-year diagonal:

```
qx(age, year) = qx_base(age) × Π over k in (2022, year] of (1 − MI(age, k, sex))
```

Projection year `t` = calendar year `2022 + t − 1`, so year 1 uses the base
table unchanged and a flat 1% scale gives `qx_base × 0.99^(t−1)` in year `t`.
Both the baseline LE (`le_pre`) and the survivor LE (`le_post`) include
improvement, so the LE change still isolates the pullforward effect. The
"SSA qx" / `orig_qx` columns in detailed output show the improved baseline.

**Importing a custom scale — CSV** (chosen over Excel: open the CSV in Excel,
edit, and "Save As → CSV"; no extra dependencies needed). Generate a
pre-filled template with `python sensitivity.py --improvement-template 1d`
(or `2d`), then run with `--improvement-table <file>`. The 1D vs 2D schema is
auto-detected from the header (any 4-digit-year column ⇒ 2D).

**1D schema** — one improvement rate per age, constant across calendar years:

```csv
age,male_improvement,female_improvement     ← or a single 'improvement' column
0,0.010,0.010                                  that applies to both sexes
1,0.010,0.010
...
```

**2D schema** — a rate per age AND calendar year (MP-scale style), wide format:

```csv
age,2023,2024,2025,...                      ← optionally an extra 'sex' column
0,0.010,0.010,0.010,...                        ('male'/'female') before the years
1,0.010,0.010,0.010,...
...
```

Rules (both schemas):
- Rates are **decimals**: `0.01` = 1% improvement/yr. Negative = deterioration.
  Any |rate| > 0.20 is rejected with an error (almost always means percents
  were entered instead of decimals).
- Ages outside the file's range use the nearest edge age; (2D) calendar years
  outside the range use the nearest edge year — years after the last column
  hold that column's rate forever, so the template's 50 year-columns
  (2023–2072) cover any projection length.
- (2D) Without a `sex` column, rates apply to both sexes; with one, each sex
  uses its own rows and a sex with no rows gets zero improvement.

---

## File Layout

```
Mortality acceleration/
├── CLAUDE.md                       # This file
├── How to run this tool.txt        # Step-by-step guide for non-technical users
├── interface.html                  # Interactive browser UI (double-click to open)
├── config.py                       # All parameters, named presets, apply_overrides()
├── main.py                         # Entry point: runs default scenario
├── sensitivity.py                  # CLI to test/sweep assumptions with no file editing
├── scenarios.py                    # Runs named scenarios; compare_all table
├── excess.py                       # CLI for the excess-calibration mode (solved pullforward)
├── improvement_template_1d.csv     # Editable 1D improvement template (age → rate)
├── improvement_template_2d.csv     # Editable 2D improvement template (age × year)
├── deploy/
│   └── index.html                  # Upload-ready copy of interface.html (refresh after edits)
├── output/                         # Created on --save; one subfolder per scenario
└── mortality_model/
    ├── __init__.py
    ├── ssa_table.py                 # SSA life table: embedded SSA 2019 + CSV loading
    ├── improvement.py               # Mortality improvement scales: flat / 1D / 2D CSV
    ├── pullforward.py               # f(t) distribution: linear, step, and exponential shapes
    ├── cohort.py                    # d(t), g(t), effective qx, mortality multiples
    ├── excess.py                    # excess-calibration: solved pullforward, calendar qx, trajectory
    └── analysis.py                  # LE calc, equivalent flat multiplier, print, CSV save
```

---

## Core Model (for reference)

### Step A — Death Distribution

`d(t)` = unconditional probability that a cohort member dies in projection year `t`. These are NOT `qx` rates; they sum to 1.0.

`d(t) = S(t-1) × qx(age + t - 1)` where `S(t) = ∏(1 - qx(age+s))` for s=0..t-1.

When mortality improvement is enabled (the default), the `qx` used here (and in
Step D's denominator) is the generational qx with improvement baked in along the
cohort's diagonal — see the Mortality Improvement section.

### Step B — Pullforward Distribution

`f(t)` = fraction of year-`t` deaths pulled into COVID. Central model assumption.

**Linear:** `f(t) = peak_fraction × max(0, 1 - (t-1) / G)`
- t=1: f = peak_fraction (maximum)
- t=G+1: f = 0

**Step:** `f(t) = peak_fraction` for `t ≤ G`, else `0` — a constant box:
"peak% of everyone who would have died within G years died during COVID
instead". Caveats: if the box spans a cohort's *entire* remaining lifetime,
removal is uniform, so with peak < 100% the survivor distribution
renormalizes back to baseline (LE change = 0 — selection needs *differential*
removal), and with peak = 100% there are no survivors and post-COVID LE is
undefined (Python raises a clear ValueError; the UI shows a warning).

**Exponential:** `f(t) = peak_fraction × exp(-decay_rate × (t-1))`

### Step C — Survivor Distribution

`d_remaining(t) = d(t) × (1 - f(t))`
`Z = Σ d_remaining(t)` = fraction of cohort surviving COVID
`g(t) = d_remaining(t) / Z` = conditional death distribution for survivors

### Step D — Mortality Multiples and Life Expectancy

`effective_qx(age+t) = g(t) / (1 - Σg(s) for s<t)`
`mortality_multiple(t) = effective_qx(age+t) / qx_SSA(age+t)` — values < 1.0 mean lower mortality for survivors

`LE = Σ S(t)` where `S(t) = 1 - Σg(s)` for s≤t

---

### Calibration mode (`excess.py` / interface "Calibrate from excess")

Runs the logic in reverse, per cohort defined by its **age in 2020** (index
`t` = calendar year 2020+t; everything per person alive at start-2020):

- **Calendar anchoring**: `qx(age, year) = qx_2019(age) × Π (1 − MI)` with the
  improvement scale applied forward from the 2019 table and *backward*
  (divided out) for years before 2019. This differs from the direct mode,
  which counts projection years from 2022.
- **Input**: cumulative excess `E` per 5-year band (fraction of one year's
  deaths at the 2020 age). Excess deaths `X = E × qx(age, 2020)`, all
  attributed to 2020.
- **Solve**: harvested deaths must repay the excess, with the deficit coming
  from the cohort *as it ages* (this mutes later years since qx rises):
  `X = Σ_{t=1..G} peak × (1 − (t−1)/G) × D_b(t)` where `D_b(t)` are baseline
  unconditional deaths. `peak` = solved share of 2021 deaths pulled into 2020.
  If `peak > 1` the case is **infeasible**: f is capped at 100% and flagged.
- **COVID path**: `D_c(t) = D_b(t) + x(t) − h(t)` inside the excess/harvest
  window (`x` = excess timing, all-2020 by default or graded over N years;
  `h(t) = f(t) × D_b(t)`), and baseline *rates* (`D_c = A_c × qx`) after the
  window — identical when feasible, graceful when capped. Multiple
  `= (D_c/A_c) / (D_b/A_b)`: spikes to 1+E in 2020, dips < 1, returns to
  exactly 1.0 after the horizon.
- **Valuation** at year V (default 2025): LE for those alive at V,
  `LE = Σ_k A(tv+k)/A(tv)`, baseline vs COVID path; the LE gain, equivalent
  flat multiplier (bisection on the improved baseline from V), and
  `alive_vs_baseline = A_c(tv)/A_b(tv)` (the not-yet-harvested deficit).
- **Trajectory** (fixed age, 2010–2035): baseline `qx(age, year)` vs the COVID
  line, where year 2020+k uses the cohort aged `age−k` in 2020 — so the
  period view stitches together many cohorts' spikes/harvests.

## Interpreting Outputs

- **Mortality multiple < 1.0**: survivors have lower mortality than the SSA table (expected — the frailest were culled)
- **LE change (years)**: additional remaining life expectancy for COVID survivors vs the SSA baseline
- **Equiv. Mort. Mult. (%)**: the flat scaling of the entire baseline qx table that reproduces the survivors' LE — e.g. 92.0 means the survivors' gain equals a uniform 8% mortality cut at every age. Lower = stronger effect.
- Effect is largest for older cohorts (more of their near-term deaths were within the pullforward window) and diminishes for younger cohorts
- The multiple returns to 1.0 at and beyond the grade-out horizon — no assumed lasting health advantage
