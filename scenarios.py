"""
Named scenarios for the mortality pullforward model.

Usage:
    python scenarios.py                     # list all scenarios
    python scenarios.py <name>              # run one scenario
    python scenarios.py compare_all         # LE-change comparison table across all scenarios
    python scenarios.py compare_all 65      # comparison table, only age 65

Available scenarios and what they test:
    default               — Baseline: linear, peak=100%, age-varying. Reference point.
    short_harvest         — Conservative: flu-harvesting view, 3yr grade-out, peak=60%
    moderate_base         — Most defensible for COVID: 7yr grade-out, peak=65%, age-varying
    elderly_concentrated  — Age-stratified: 85+ gets steep (3yr), <65 gets flat (15yr)
    long_harvest          — Aggressive upper bound: 15yr grade-out, peak=50%
    exponential_moderate  — Fat-tail, no hard cutoff, peak=70%, decay=0.4
    exponential_long_tail — Fat-tail, slow decay, peak=80%, decay=0.2
    sickest_only          — Very short (2yr), high peak (90%): only the imminently dying
    uniform_benchmark     — Flat 10yr grade-out, peak=50%, no age-variation
"""
from __future__ import annotations
import sys

from config import (
    describe_improvement,
    default_config,
    short_harvest_config,
    moderate_base_config,
    elderly_concentrated_config,
    long_harvest_config,
    exponential_moderate_config,
    exponential_long_tail_config,
    sickest_only_config,
    uniform_benchmark_config,
)
from mortality_model.analysis import run_all_cohorts, print_results, print_summary_table

import pandas as pd

SCENARIOS = {
    "default":               (default_config,               "Baseline: linear, peak=100%, age-varying grade-out"),
    "short_harvest":         (short_harvest_config,         "Conservative: 3yr grade-out, peak=60% -- flu harvesting view"),
    "moderate_base":         (moderate_base_config,         "Moderate: 7yr grade-out, peak=65%, age-varying -- most defensible"),
    "elderly_concentrated":  (elderly_concentrated_config,  "Age-stratified: steep for 85+ (3yr), flat for <65 (15yr), peak=70%"),
    "long_harvest":          (long_harvest_config,          "Aggressive: 15yr grade-out, peak=50% -- upper bound"),
    "exponential_moderate":  (exponential_moderate_config,  "Exponential fat-tail, decay=0.4, peak=70%"),
    "exponential_long_tail": (exponential_long_tail_config, "Exponential slow decay, decay=0.2, peak=80%"),
    "sickest_only":          (sickest_only_config,          "Only imminently dying: 2yr grade-out, peak=90%"),
    "uniform_benchmark":     (uniform_benchmark_config,     "Clean benchmark: 10yr, peak=50%, no age-variation"),
}


def run_scenario(name: str) -> list:
    if name not in SCENARIOS:
        print(f"Unknown scenario: {name!r}")
        print(f"Available: {list(SCENARIOS.keys())}")
        sys.exit(1)

    factory, description = SCENARIOS[name]
    cfg = factory()
    pf = cfg.pullforward

    print(f"\n{'='*65}")
    print(f"Scenario: {name}")
    print(f"  {description}")
    print(f"  Shape          : {pf.shape}")
    print(f"  Peak fraction  : {pf.peak_fraction:.0%}")
    print(f"  Age-varying    : {pf.age_varying}")
    if pf.shape in ("linear", "step"):
        if pf.age_varying:
            print(f"  Grade-out bands: {pf.age_grade_out}")
        else:
            print(f"  Grade-out years: {pf.default_grade_out_years}")
    else:
        print(f"  Decay rate     : {pf.exponential_decay_rate}")
    print(f"  Improvement    : {describe_improvement(cfg)}")
    print(f"{'='*65}")

    results = run_all_cohorts(cfg)
    save = "--save" in sys.argv
    print_results(results, scenario_name=name, save=save)
    return results


def compare_all(filter_age: int = None) -> None:
    """
    Run all scenarios and print LE-change side by side.
    Optionally filter to a single age (e.g. compare_all(65)).
    """
    all_rows = []
    for scenario_name, (factory, _) in SCENARIOS.items():
        cfg = factory()
        results = run_all_cohorts(cfg)
        for r in results:
            if filter_age is not None and r["starting_age"] != filter_age:
                continue
            all_rows.append({
                "Scenario":       scenario_name,
                "Sex":            r["sex"].capitalize(),
                "Age":            r["starting_age"],
                "LE Pre (yrs)":   round(r["le_pre"], 2),
                "LE Post (yrs)":  round(r["le_post"], 2),
                "LE Change (yrs)": round(r["le_change"], 2),
                "LE Change (%)":  round(r["le_pct_change"], 1),
                "Equiv Mult (%)": round(r["equiv_flat_multiplier"] * 100, 1),
            })

    df = pd.DataFrame(all_rows)

    print("\n=== Scenario Comparison: LE Change (years) by Age/Sex ===")
    pivot = df.pivot_table(
        index=["Sex", "Age"],
        columns="Scenario",
        values="LE Change (yrs)",
    )[list(SCENARIOS.keys())]  # preserve scenario order
    print(pivot.to_string())

    print("\n=== Scenario Comparison: LE Change (%) by Age/Sex ===")
    pivot_pct = df.pivot_table(
        index=["Sex", "Age"],
        columns="Scenario",
        values="LE Change (%)",
    )[list(SCENARIOS.keys())]
    print(pivot_pct.to_string())

    print("\n=== Scenario Comparison: Equivalent Flat Mortality Multiplier (% of baseline qx) ===")
    print("(flat scaling of the entire qx table that reproduces the survivors' LE)")
    pivot_mult = df.pivot_table(
        index=["Sex", "Age"],
        columns="Scenario",
        values="Equiv Mult (%)",
    )[list(SCENARIOS.keys())]
    print(pivot_mult.to_string())


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Available scenarios:\n")
        for name, (_, desc) in SCENARIOS.items():
            print(f"  {name:30s} {desc}")
        print("\nUsage:")
        print("  python scenarios.py <scenario_name>")
        print("  python scenarios.py compare_all")
        print("  python scenarios.py compare_all 65   # filter to age 65")
        sys.exit(0)

    arg = sys.argv[1]
    if arg == "compare_all":
        filter_age = int(sys.argv[2]) if len(sys.argv) > 2 else None
        compare_all(filter_age)
    else:
        run_scenario(arg)
