# Mortality Pullforward Model — User Guide

This model estimates how much COVID-19 increased the remaining life expectancy of survivors, by pulling the deaths of the sickest people (who would have died soonest anyway) into the COVID period. The result: the post-COVID survivor population is healthier on average than the pre-COVID population at the same age.

---

## Setup (one-time)

You need Python installed. Then, in a terminal in this folder, install the required packages:

```
pip install -r requirements.txt
```

---

## Running the Model

Open a terminal, navigate to this folder, and use one of these commands:

### Run the default scenario and print results
```
python main.py
```

### Run the default scenario and save results to CSV files
```
python main.py --save
```

### Run a specific named scenario
```
python scenarios.py moderate_base
python scenarios.py moderate_base --save
```

### See all available scenario names
```
python scenarios.py
```

### Compare all scenarios side by side
```
python scenarios.py compare_all
python scenarios.py compare_all 65        (filter to age 65 only)
```

---

## Where Results Are Saved

When you run with `--save`, results appear in the `output/` folder:

```
output/
└── moderate_base/
    ├── summary.csv                        ← main results: LE before/after, change, % survived
    ├── multiples_male_age65.csv           ← year-by-year detail for 65-year-old males
    ├── multiples_female_age65.csv
    └── ...                                ← one detail file per age and sex
```

**summary.csv columns:**

| Column | Meaning |
|---|---|
| `sex` | male or female |
| `age` | starting age of the cohort |
| `le_pre_covid` | remaining life expectancy before COVID (from SSA table), in years |
| `le_post_covid` | remaining life expectancy for COVID survivors, in years |
| `le_change_years` | how many additional years of LE survivors gained |
| `le_change_pct` | same, as a percentage |
| `survivor_fraction` | fraction of the original cohort that survived COVID |

**multiples_[sex]_age[N].csv columns:**

| Column | Meaning |
|---|---|
| `age` | actual age in that projection year |
| `years_from_now` | years after the end of COVID (2022) |
| `orig_qx` | annual death probability from the SSA table |
| `pullforward_fraction` | f(t): fraction of that year's deaths pulled into COVID |
| `d_original` | unconditional probability of dying in this year (pre-COVID) |
| `g_survivor` | unconditional probability of dying in this year (for COVID survivors) |
| `effective_qx` | implied annual death probability for survivors |
| `mortality_multiple` | effective_qx ÷ orig_qx. Values below 1.0 mean survivors are healthier than average |

---

## How to Change the Assumptions

All inputs are controlled in **`config.py`**. Open that file in any text editor.

### The most important parameters

**`peak_fraction`** — How severe is the pullforward at its peak?

At `t=1` (people who would have died in the first year after COVID), what fraction of those deaths were actually pulled into COVID? This is the single biggest lever.

- `1.0` = 100% of imminent deaths were pulled in (aggressive, used as theoretical ceiling)
- `0.65` = 65% (moderate, most defensible)
- `0.30` = 30% (conservative)

**`default_grade_out_years`** (for linear shape) — How many years out does the pullforward reach?

- `3` = COVID only pulled forward deaths that were within ~3 years (conservative, like flu harvesting)
- `7`–`10` = moderate
- `15` = aggressive upper bound

**`shape`** — Shape of the pullforward curve:

- `"linear"` — grades from peak down to zero over `grade_out_years`, then stops completely
- `"exponential"` — never fully stops, but decays quickly (fat tail). Controlled by `exponential_decay_rate`

**`age_varying`** — Should older cohorts have a steeper pullforward?

- `True` = yes: an 85-year-old's pullforward grades out in 3–5 years, a 45-year-old's grades out in 10+ years. This reflects that COVID was most impactful among the frailest elderly.
- `False` = same grade-out for all ages

### Example: creating a custom scenario

In `config.py`, copy any existing preset function and change the parameters:

```python
def my_custom_config() -> ModelConfig:
    pf = PullforwardConfig(
        shape="linear",
        peak_fraction=0.55,          # 55% of year-1 deaths pulled in
        default_grade_out_years=6,   # effect grades out over 6 years
        age_varying=True,
        age_grade_out={
            (0, 60):   10,
            (60, 70):   6,
            (70, 80):   4,
            (80, 200):  3,
        },
    )
    return ModelConfig(pullforward=pf)
```

Then in `scenarios.py`, add it to the `SCENARIOS` dictionary near the top of the file:

```python
SCENARIOS = {
    ...
    "my_custom": (my_custom_config, "My custom assumptions"),
}
```

Then run it:
```
python scenarios.py my_custom --save
```

---

## Available Named Scenarios

These are pre-built to represent a range of views on how severe the pullforward was:

| Name | Peak | Grade-out | Notes |
|---|---|---|---|
| `default` | 100% | 10 yrs, age-varying | Theoretical upper bound on LE gain |
| `short_harvest` | 60% | 3 yrs, uniform | Conservative; similar to flu harvesting research |
| `moderate_base` | 65% | 7 yrs, age-varying | Most defensible for COVID specifically |
| `elderly_concentrated` | 70% | 3–15 yrs, strong age-variation | Concentrated impact on 85+; mild for <65 |
| `long_harvest` | 50% | 15 yrs, age-varying | Aggressive upper bound |
| `exponential_moderate` | 70% | Fat tail, decay=0.4 | No hard cutoff; moderate decay |
| `exponential_long_tail` | 80% | Fat tail, decay=0.2 | Slow decay; broad long-run effect |
| `sickest_only` | 90% | 2 yrs, uniform | Only the imminently dying; nearly all of them |
| `uniform_benchmark` | 50% | 10 yrs, no age-variation | Clean benchmark |

**Recommended starting point:** `moderate_base` and `short_harvest` bracket the most defensible range.

---

## How to Use Your Own Mortality Table

The model uses a built-in approximation of the SSA 2019 Period Life Table by default. To use a different table:

1. Prepare a CSV file with exactly these columns: `age`, `male_qx`, `female_qx`
   - One row per age, from 0 to 119
   - `qx` values are the annual probability of death (between 0 and 1)
   - Age 119 must have `qx = 1.0`

2. In `config.py`, set the path in `ModelConfig`:
   ```python
   ssa_table_path: Optional[str] = "data/my_table.csv"
   ```

---

## How to Interpret the Results

**Life expectancy change:** A positive number means COVID survivors at that age have a longer remaining life expectancy than the SSA table would predict. For example, `le_change_years = 1.2` for a 65-year-old male means survivors are expected to live 1.2 years longer than the baseline table suggests.

**Mortality multiple:** The ratio of survivors' mortality to the SSA table mortality at each future age. A value of `0.70` means survivors have 30% lower mortality than the table at that age. The multiple starts below 1.0 (survivors are healthier) and returns to 1.0 at and beyond the grade-out horizon (no assumed lasting health advantage after that point).

**Survivor fraction:** The fraction of the original cohort that made it through COVID. A value of `0.89` for 65-year-old males means 89% survived, and the 11% who died were disproportionately the sicker members of the cohort.

**General pattern to expect:**
- Older cohorts show bigger LE gains in percentage terms, because more of their near-term deaths fell within the pullforward window
- Males show larger absolute LE gains than females at the same age, because male baseline mortality is higher
- The effect shrinks to zero beyond the grade-out horizon by construction
