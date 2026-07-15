"""
Mortality Pullforward Model — Main Entry Point

Runs the default scenario and prints:
  1. Summary table: pre/post-COVID life expectancy and change by age/sex
  2. Mortality multiples vs SSA table at key future ages

Edit config.py to change assumptions, or run scenarios.py for named presets.
"""
import sys
from config import default_config, describe_improvement
from mortality_model.analysis import run_all_cohorts, print_results, run_cohort_analysis
from mortality_model.ssa_table import load_ssa_table
import numpy as np
import pandas as pd


def main():
    cfg = default_config()

    print(f"Mortality Pullforward Model")
    print(f"  Reference year (end of COVID): {cfg.covid_end_year}")
    print(f"  Pullforward shape: {cfg.pullforward.shape}")
    print(f"  Age-varying grade-out: {cfg.pullforward.age_varying}")
    if cfg.pullforward.age_varying:
        print(f"  Grade-out by age band: {cfg.pullforward.age_grade_out}")
    else:
        print(f"  Default grade-out: {cfg.pullforward.default_grade_out_years} years")
    print(f"  Mortality improvement: {describe_improvement(cfg)}")

    save = "--save" in sys.argv
    results = run_all_cohorts(cfg)
    print_results(results, scenario_name="default", save=save)

    # --- Optional: detailed per-cohort output for one specific age/sex ---
    print_detailed_cohort(results, target_age=65, target_sex="male")
    print_detailed_cohort(results, target_age=65, target_sex="female")


def print_detailed_cohort(results: list, target_age: int, target_sex: str) -> None:
    """Print year-by-year death distribution and mortality multiples for one cohort."""
    match = [r for r in results if r["starting_age"] == target_age and r["sex"] == target_sex]
    if not match:
        return
    r = match[0]
    df = r["multiples_df"]

    print(f"\n--- Detailed: {target_sex.capitalize()}, Age {target_age} ---")
    print(f"  Survivor fraction (% survived COVID): {r['survivor_fraction']*100:.1f}%")
    print(f"  Pre-COVID LE:  {r['le_pre']:.2f} years")
    print(f"  Post-COVID LE: {r['le_post']:.2f} years")
    print(f"  LE Change:     {r['le_change']:+.2f} years  ({r['le_pct_change']:+.1f}%)")

    # Show first 25 years in detail
    view = df[df["years_from_now"] <= 25].copy()
    view = view[["years_from_now", "age", "orig_qx", "pullforward_fraction",
                 "d_original", "g_survivor", "effective_qx", "mortality_multiple"]].copy()
    view.columns = ["Year", "Age", "SSA qx", "Pullforward f(t)",
                    "d(t) Original", "g(t) Survivor", "Eff. qx", "Mult."]
    view = view.set_index("Year")

    # Format floats
    fmt = {
        "SSA qx": "{:.5f}",
        "Pullforward f(t)": "{:.3f}",
        "d(t) Original": "{:.5f}",
        "g(t) Survivor": "{:.5f}",
        "Eff. qx": "{:.5f}",
        "Mult.": "{:.4f}",
    }
    print(view.to_string(float_format=lambda x: f"{x:.5f}"))


if __name__ == "__main__":
    main()
