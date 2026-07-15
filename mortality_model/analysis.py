"""
Life expectancy calculations and result formatting.

Runs the full cohort model for each (age, sex) pair and produces:
  - Pre-COVID life expectancy (from SSA table alone)
  - Post-COVID life expectancy (for COVID survivors)
  - Life expectancy change (years gained)
  - Mortality multiples by age
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from config import ModelConfig
from mortality_model.ssa_table import load_ssa_table, get_qx
from mortality_model.improvement import (
    ImprovementScale,
    load_improvement_scale,
    build_cohort_table,
)
from mortality_model.cohort import (
    compute_death_distribution,
    apply_pullforward,
    compute_survivor_distribution,
    compute_effective_qx,
    compute_mortality_multiples,
)


def compute_life_expectancy(death_dist: np.ndarray) -> float:
    """
    Compute remaining life expectancy from a death distribution.

    Uses: E[T] = sum_{t=0}^{inf} S(t)  where S(t) = P(T > t).

    death_dist[i] = P(die in year i+1), 0-indexed.
    Returns expected remaining years.
    """
    n = len(death_dist)
    # S(0) = 1 (alive at start), S(i) = S(i-1) - d(i-1)
    survival = np.ones(n + 1)
    for i in range(1, n + 1):
        survival[i] = survival[i - 1] - death_dist[i - 1]

    # E[T] = S(0) + S(1) + ... + S(n-1)   (S(n) is essentially 0 by construction)
    return float(survival[:-1].sum())


def _le_from_qx(qx: np.ndarray) -> float:
    """Remaining LE implied by a vector of annual qx (terminal qx forced to 1)."""
    q = np.clip(qx, 0.0, 1.0).copy()
    q[-1] = 1.0
    surv = np.cumprod(1.0 - q)
    d = np.empty_like(q)
    d[0] = q[0]
    d[1:] = surv[:-1] * q[1:]
    total = d.sum()
    if total > 0:
        d /= total
    return compute_life_expectancy(d)


def solve_equivalent_flat_multiplier(
    target_le: float,
    starting_age: int,
    sex: str,
    table: pd.DataFrame,
    max_age: int = 119,
) -> float:
    """
    Find the single multiplier m such that scaling qx at EVERY age of the
    baseline table by m reproduces target_le for this cohort.

    This expresses the survivors' LE gain as an equivalent flat percentage
    adjustment to the whole mortality table: m = 0.90 means the survivors'
    life expectancy is what you would get by cutting every qx by 10%.

    LE is strictly decreasing in m, so a bisection converges cleanly.
    """
    qx = np.array([get_qx(table, a, sex) for a in range(starting_age, max_age + 1)])

    lo, hi = 1e-6, 1.0
    # Survivors normally live longer than baseline (m < 1), but allow m > 1
    # in case a scale/assumption combination produces the reverse.
    while _le_from_qx(qx * hi) > target_le and hi < 64.0:
        lo = hi
        hi *= 2.0
    if _le_from_qx(qx * lo) < target_le:
        return np.nan  # target above even the near-zero-mortality LE

    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if _le_from_qx(qx * mid) > target_le:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def run_cohort_analysis(
    starting_age: int,
    sex: str,
    cfg: ModelConfig,
    ssa_table: pd.DataFrame,
    improvement_scale: ImprovementScale = None,
) -> dict:
    """
    Full analysis for one (age, sex) cohort.

    Returns a dict with:
      - starting_age, sex
      - survivor_fraction: fraction of cohort that survived COVID
      - le_pre: pre-COVID remaining LE (years)
      - le_post: post-COVID remaining LE for survivors (years)
      - le_change: le_post - le_pre
      - le_pct_change: percentage change
      - equiv_flat_multiplier: flat qx multiplier on the whole baseline table
        that reproduces le_post (< 1.0 means survivors act like a table with
        uniformly lower mortality)
      - multiples_df: DataFrame of mortality multiples by future age

    If improvement_scale is given, mortality improvement is baked into the
    baseline table along the cohort's age/calendar-year diagonal first, so
    both le_pre and le_post (and the multiples' denominator) include it.
    """
    # Step 0: Bake mortality improvement into this cohort's baseline table
    if improvement_scale is not None:
        ssa_table = build_cohort_table(
            ssa_table, starting_age, improvement_scale,
            cfg.covid_end_year, cfg.max_age,
        )

    # Step 1: Original death distribution
    d = compute_death_distribution(starting_age, sex, ssa_table, cfg.max_age)

    # Step 2: Apply pullforward
    d_remaining, f = apply_pullforward(d, starting_age, sex, cfg)

    # Step 3: Survivor conditional distribution
    g, survivor_fraction = compute_survivor_distribution(d_remaining)

    # Step 4: Effective qx and mortality multiples
    eff_qx, ages = compute_effective_qx(g, starting_age, ssa_table)
    multiples = compute_mortality_multiples(eff_qx, starting_age, sex, ssa_table)

    # Step 5: Life expectancy
    le_pre = compute_life_expectancy(d)
    le_post = compute_life_expectancy(g)
    le_change = le_post - le_pre

    # Step 6: Equivalent flat mortality multiplier — the uniform qx scaling of
    # the whole baseline table that would deliver the same LE as the survivors
    equiv_mult = solve_equivalent_flat_multiplier(
        le_post, starting_age, sex, ssa_table, cfg.max_age
    )

    # Build per-year DataFrame
    n = len(d)
    orig_qx_series = np.array([
        float(ssa_table.loc[min(starting_age + i, cfg.max_age), f"{sex}_qx"])
        for i in range(n)
    ])
    multiples_df = pd.DataFrame({
        "age": ages,
        "years_from_now": np.arange(1, n + 1),
        "orig_qx": orig_qx_series,
        "pullforward_fraction": f,
        "d_original": d,
        "d_remaining": d_remaining,
        "g_survivor": g,
        "effective_qx": eff_qx,
        "mortality_multiple": multiples,
    })

    return {
        "starting_age": starting_age,
        "sex": sex,
        "survivor_fraction": survivor_fraction,
        "le_pre": le_pre,
        "le_post": le_post,
        "le_change": le_change,
        "le_pct_change": le_change / le_pre * 100 if le_pre > 0 else np.nan,
        "equiv_flat_multiplier": equiv_mult,
        "multiples_df": multiples_df,
    }


def run_all_cohorts(cfg: ModelConfig) -> list[dict]:
    """Run analysis for all (age, sex) pairs in cfg.analysis_ages / cfg.analysis_sexes."""
    ssa_table = load_ssa_table(cfg.ssa_table_path, cfg.max_age)
    scale = None
    if cfg.improvement_enabled:
        scale = load_improvement_scale(cfg.improvement_table_path, cfg.improvement_rate)
    results = []
    for sex in cfg.analysis_sexes:
        for age in cfg.analysis_ages:
            result = run_cohort_analysis(age, sex, cfg, ssa_table, scale)
            results.append(result)
    return results


# ---------------------------------------------------------------------------
# Formatted output helpers
# ---------------------------------------------------------------------------

def print_summary_table(results: list[dict]) -> None:
    """Print a high-level summary: LE before/after and change."""
    rows = []
    for r in results:
        rows.append({
            "Sex": r["sex"].capitalize(),
            "Age": r["starting_age"],
            "LE Pre-COVID (yrs)": round(r["le_pre"], 2),
            "LE Post-COVID (yrs)": round(r["le_post"], 2),
            "LE Change (yrs)": round(r["le_change"], 2),
            "LE Change (%)": round(r["le_pct_change"], 1),
            "Equiv. Mort. Mult. (%)": round(r["equiv_flat_multiplier"] * 100, 1),
            "% Survived COVID": round(r["survivor_fraction"] * 100, 1),
        })
    df = pd.DataFrame(rows).set_index(["Sex", "Age"])
    print("\n=== Life Expectancy Impact of COVID Pullforward ===")
    print(df.to_string())
    print("\n  Equiv. Mort. Mult. (%): flat scaling of the ENTIRE baseline qx table")
    print("  that reproduces the survivors' LE -- e.g. 92.0 means survivors gain the")
    print("  same LE as a uniform 8% cut in mortality at every age.")


def print_mortality_multiples(results: list[dict], sample_ages_offset: list[int] = None) -> None:
    """
    Print mortality multiples at selected future ages for each cohort.

    sample_ages_offset: years into the future to show (default: [0, 5, 10, 15, 20])
    """
    if sample_ages_offset is None:
        sample_ages_offset = [0, 5, 10, 15, 20]

    print("\n=== Mortality Multiples vs SSA Table (post-COVID survivors) ===")
    print("(Values < 1.0 mean survivors have LOWER mortality than SSA table)\n")

    for r in results:
        df = r["multiples_df"]
        label = f"{r['sex'].capitalize()}, starting age {r['starting_age']}"
        print(f"  {label}:")
        header = f"  {'Years From Now':>16}  {'Age':>5}  {'SSA qx':>9}  {'Eff. qx':>9}  {'Multiple':>9}"
        print(header)
        for offset in sample_ages_offset:
            row = df[df["years_from_now"] == offset + 1]
            if row.empty:
                continue
            row = row.iloc[0]
            mult_str = f"{row['mortality_multiple']:.3f}" if not np.isnan(row["mortality_multiple"]) else "  N/A"
            eff_str = f"{row['effective_qx']:.4f}" if not np.isnan(row["effective_qx"]) else "  N/A"
            print(
                f"  {offset:>16}  {int(row['age']):>5}  "
                f"{row['orig_qx']:>9.4f}  {eff_str:>9}  {mult_str:>9}"
            )
        print()


def save_results(results: list[dict], scenario_name: str, output_dir: str = "output") -> None:
    """
    Save results to CSVs in output/<scenario_name>/.

    Files written:
      summary.csv         — one row per (age, sex): LE pre/post/change
      multiples_<sex>_age<N>.csv — year-by-year mortality multiples per cohort
    """
    import os
    folder = os.path.join(output_dir, scenario_name)
    os.makedirs(folder, exist_ok=True)

    # Summary table
    rows = []
    for r in results:
        rows.append({
            "sex": r["sex"],
            "age": r["starting_age"],
            "le_pre_covid": round(r["le_pre"], 4),
            "le_post_covid": round(r["le_post"], 4),
            "le_change_years": round(r["le_change"], 4),
            "le_change_pct": round(r["le_pct_change"], 2),
            "equiv_flat_qx_multiplier": round(r["equiv_flat_multiplier"], 6),
            "survivor_fraction": round(r["survivor_fraction"], 4),
        })
    summary_df = pd.DataFrame(rows)
    summary_path = os.path.join(folder, "summary.csv")
    summary_df.to_csv(summary_path, index=False)
    print(f"  Saved: {summary_path}")

    # Per-cohort detail
    for r in results:
        fname = f"multiples_{r['sex']}_age{r['starting_age']}.csv"
        fpath = os.path.join(folder, fname)
        r["multiples_df"].to_csv(fpath, index=False)

    print(f"  Saved: {folder}/multiples_*.csv  ({len(results)} files)")


def print_results(results: list[dict], scenario_name: str = None, save: bool = False) -> None:
    """Print full formatted output, and optionally save to output/."""
    print_summary_table(results)
    print_mortality_multiples(results)
    if save and scenario_name:
        save_results(results, scenario_name)
