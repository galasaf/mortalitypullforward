"""
sensitivity.py - Test how the COVID pullforward results respond to changing
assumptions, WITHOUT editing any files.

Run `python sensitivity.py` (no arguments) to see the worked examples below.

----------------------------------------------------------------------------
MODE 1 -- SINGLE RUN: try one combination of assumptions
----------------------------------------------------------------------------
  python sensitivity.py --peak 0.6 --grade-out 7
  python sensitivity.py --shape exponential --decay 0.3 --age 65
  python sensitivity.py --base moderate_base --peak 0.8        # tweak a preset

----------------------------------------------------------------------------
MODE 2 -- SWEEP: vary ONE knob across a range, see how LE change responds
----------------------------------------------------------------------------
  python sensitivity.py --sweep peak       --range 0.3 0.8 0.1
  python sensitivity.py --sweep grade-out  --values 2 3 5 7 10 15
  python sensitivity.py --sweep decay --range 0.1 0.5 0.1 --shape exponential

----------------------------------------------------------------------------
MODE 3 -- MORTALITY IMPROVEMENT (default: flat 1%/yr at every age and year)
----------------------------------------------------------------------------
  python sensitivity.py --improvement-rate 0.015          # flat 1.5%/yr
  python sensitivity.py --no-improvement                  # static life table
  python sensitivity.py --improvement-table my_scale.csv  # 1D or 2D CSV
  python sensitivity.py --improvement-template 1d         # write a template & exit
  python sensitivity.py --improvement-template 2d
  python sensitivity.py --sweep improvement --range 0.0 0.02 0.005 --age 65

----------------------------------------------------------------------------
COMMON OPTIONS (work in both modes)
----------------------------------------------------------------------------
  --age 65              focus on one starting age (default: all configured ages)
  --sex male|female     focus on one sex (default: both)
  --base <scenario>     start from a named preset instead of model defaults
  --save                also write CSVs to output/<name>/   (single-run only)

Knobs you can override / sweep:
  peak        peak_fraction   (0..1)  fraction of year-1 deaths pulled into COVID
  grade-out   horizon in years (linear: ramps to zero there; step: constant
              full effect inside it, zero after; turns off age bands)
  decay       exponential decay rate  (only matters when --shape exponential)
  shape       linear | step | exponential
  improvement flat annual mortality improvement rate (0.01 = 1%/yr)
"""
from __future__ import annotations

import argparse
import sys

import pandas as pd

from config import (
    ModelConfig,
    apply_overrides,
    describe_config,
)
from mortality_model.analysis import run_all_cohorts, print_results

# Reuse the named presets as optional starting points for --base.
from scenarios import SCENARIOS

# Map the friendly sweep name to the apply_overrides keyword.
SWEEPABLE = {
    "peak": "peak",
    "grade-out": "grade_out",
    "grade_out": "grade_out",
    "decay": "decay",
    "improvement": "improvement_rate",
}


def base_config(name: str | None) -> ModelConfig:
    """Return a fresh ModelConfig, optionally seeded from a named preset."""
    if name is None:
        return ModelConfig()
    if name not in SCENARIOS:
        print(f"Unknown --base scenario: {name!r}")
        print(f"Available: {', '.join(SCENARIOS)}")
        sys.exit(1)
    factory, _ = SCENARIOS[name]
    return factory()


def resolve_scope(args):
    """Translate --age / --sex into analysis_ages / analysis_sexes overrides."""
    ages = [args.age] if args.age is not None else None
    sexes = [args.sex] if args.sex is not None else None
    return ages, sexes


def cmd_single(args):
    """One run with the given overrides; prints the standard report."""
    ages, sexes = resolve_scope(args)
    cfg = base_config(args.base)
    apply_overrides(
        cfg,
        peak=args.peak,
        grade_out=args.grade_out,
        shape=args.shape,
        decay=args.decay,
        ages=ages,
        sexes=sexes,
        improvement_rate=args.improvement_rate,
        improvement_table=args.improvement_table,
        no_improvement=args.no_improvement,
    )

    print("=" * 65)
    print("Single run" + (f" (base: {args.base})" if args.base else ""))
    print(describe_config(cfg))
    print("=" * 65)

    try:
        results = run_all_cohorts(cfg)
    except ValueError as e:
        print(f"\nERROR: {e}")
        sys.exit(1)
    name = args.base or "adhoc"
    print_results(results, scenario_name=name, save=args.save)


def parse_sweep_values(args) -> list[float]:
    """Build the list of values to sweep from --values or --range LOW HIGH STEP."""
    if args.values is not None:
        return list(args.values)
    if args.range is not None:
        low, high, step = args.range
        if step <= 0:
            print("--range STEP must be positive.")
            sys.exit(1)
        vals, v = [], low
        # Use a small epsilon so the high endpoint is included despite FP drift.
        while v <= high + step * 1e-9:
            vals.append(round(v, 10))
            v += step
        return vals
    print("Sweep needs either --values v1 v2 ... or --range LOW HIGH STEP")
    sys.exit(1)


def cmd_sweep(args):
    """Vary one knob across values; print a tidy LE-change table."""
    knob = SWEEPABLE.get(args.sweep)
    if knob is None:
        print(f"Cannot sweep {args.sweep!r}. Sweepable: {', '.join(sorted(set(SWEEPABLE)))}")
        sys.exit(1)

    values = parse_sweep_values(args)
    if knob == "grade_out":
        values = [int(round(v)) for v in values]  # grade-out is whole years
    ages, sexes = resolve_scope(args)

    # grade-out sweeps imply a linear, flat (non-age-varying) curve.
    forced_shape = args.shape
    if knob == "grade_out" and forced_shape is None:
        forced_shape = "linear"
    if knob == "decay" and forced_shape is None:
        forced_shape = "exponential"

    print("=" * 65)
    print(f"SWEEP of '{args.sweep}' over {values}")
    if args.base:
        print(f"  base preset    : {args.base}")
    print(f"  shape          : {forced_shape or '(from base/default)'}")
    if args.age is not None:
        print(f"  age            : {args.age}")
    if args.sex is not None:
        print(f"  sex            : {args.sex}")
    print("=" * 65)

    rows = []
    for val in values:
        cfg = base_config(args.base)
        overrides = dict(shape=forced_shape, ages=ages, sexes=sexes)
        overrides[knob] = val
        # Keep any non-swept knobs the user also passed on the command line.
        if knob != "peak" and args.peak is not None:
            overrides["peak"] = args.peak
        if knob != "decay" and args.decay is not None:
            overrides["decay"] = args.decay
        if knob != "grade_out" and args.grade_out is not None:
            overrides["grade_out"] = args.grade_out
        if knob != "improvement_rate":
            overrides["improvement_rate"] = args.improvement_rate
            overrides["improvement_table"] = args.improvement_table
            overrides["no_improvement"] = args.no_improvement
        apply_overrides(cfg, **overrides)

        try:
            cohort_results = run_all_cohorts(cfg)
        except ValueError as e:
            print(f"  [skipping {args.sweep}={val}: {e}]")
            continue
        for r in cohort_results:
            rows.append({
                args.sweep: val,
                "Sex": r["sex"].capitalize(),
                "Age": r["starting_age"],
                "LE Pre": round(r["le_pre"], 2),
                "LE Post": round(r["le_post"], 2),
                "LE Change (yrs)": round(r["le_change"], 2),
                "LE Change (%)": round(r["le_pct_change"], 1),
                "Equiv Mult (%)": round(r["equiv_flat_multiplier"] * 100, 1),
                "% Survived": round(r["survivor_fraction"] * 100, 1),
            })

    df = pd.DataFrame(rows)

    # If a single age+sex, a flat table reads best. Otherwise pivot so each
    # swept value is a column and each (Sex, Age) is a row.
    n_ages = df["Age"].nunique()
    n_sex = df["Sex"].nunique()
    if n_ages == 1 and n_sex == 1:
        print("\n" + df.drop(columns=["Sex", "Age"]).to_string(index=False))
    else:
        print("\n=== LE Change (years) - rows: Sex/Age, columns: "
              f"{args.sweep} ===")
        pivot = df.pivot_table(index=["Sex", "Age"], columns=args.sweep,
                               values="LE Change (yrs)")
        print(pivot.to_string())
        print("\n=== LE Change (%) ===")
        pivot_pct = df.pivot_table(index=["Sex", "Age"], columns=args.sweep,
                                   values="LE Change (%)")
        print(pivot_pct.to_string())
        print("\n=== Equivalent flat mortality multiplier (% of baseline qx) ===")
        pivot_mult = df.pivot_table(index=["Sex", "Age"], columns=args.sweep,
                                    values="Equiv Mult (%)")
        print(pivot_mult.to_string())


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Test sensitivity of the mortality pullforward model "
                    "without editing files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    # Override knobs (used in both modes)
    p.add_argument("--peak", type=float, help="peak_fraction, 0..1")
    p.add_argument("--grade-out", type=int, dest="grade_out",
                   help="linear grade-out years (turns off age bands)")
    p.add_argument("--shape", choices=["linear", "step", "exponential"])
    p.add_argument("--decay", type=float, help="exponential decay rate")
    p.add_argument("--base", help="start from a named scenario preset")
    p.add_argument("--age", type=int, help="focus on one starting age")
    p.add_argument("--sex", choices=["male", "female"], help="focus on one sex")
    p.add_argument("--save", action="store_true",
                   help="write CSVs (single-run mode only)")

    # Mortality improvement
    p.add_argument("--improvement-rate", type=float, dest="improvement_rate",
                   help="flat annual mortality improvement rate (0.01 = 1%%/yr; default 0.01)")
    p.add_argument("--improvement-table", dest="improvement_table",
                   help="path to a 1D or 2D improvement CSV (schema auto-detected)")
    p.add_argument("--no-improvement", action="store_true", dest="no_improvement",
                   help="disable mortality improvement (static life table)")
    p.add_argument("--improvement-template", choices=["1d", "2d"],
                   dest="improvement_template",
                   help="write improvement_template_<1d|2d>.csv (pre-filled with 1%%) and exit")

    # Sweep mode
    p.add_argument("--sweep", help="knob to vary: peak | grade-out | decay | improvement")
    p.add_argument("--values", type=float, nargs="+",
                   help="explicit values to sweep, e.g. --values 2 5 10")
    p.add_argument("--range", type=float, nargs=3, metavar=("LOW", "HIGH", "STEP"),
                   help="sweep LOW..HIGH inclusive, in steps of STEP")
    return p


def main():
    if len(sys.argv) == 1:
        print(__doc__)
        return
    args = build_parser().parse_args()
    if args.improvement_template:
        from mortality_model.improvement import write_template_1d, write_template_2d
        writer = write_template_1d if args.improvement_template == "1d" else write_template_2d
        path = writer()
        print(f"Wrote {path} (pre-filled with 1% at every age/year).")
        print("Edit the rates (decimals: 0.01 = 1%/yr), then run:")
        print(f"  python sensitivity.py --improvement-table {path}")
        return
    if args.sweep:
        cmd_sweep(args)
    else:
        cmd_single(args)


if __name__ == "__main__":
    main()
