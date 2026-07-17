"""
excess.py - the unified COVID pullforward model (calendar-anchored).

The pullforward has two equivalent descriptions tied together by one
conservation equation (excess deaths in 2020 = deaths harvested from later
years), so you can drive the model from EITHER end:

  * Assume the EXCESS (default): the cumulative excess mortality per age
    group (% of one normal year's deaths, in 5-year bands) -> the model
    SOLVES the pullforward peak (the share of 2021's deaths pulled into
    2020).
  * Assume the PULLFORWARD (--peak): the peak plus a grade-away shape ->
    the model reports the cumulative excess this IMPLIES.

Standing at the valuation year (default 2025), it reports the remaining
mortality dip and the life-expectancy gain for people still alive.

Run `python excess.py` (no arguments) for these worked examples.

----------------------------------------------------------------------------
BASIC RUNS
----------------------------------------------------------------------------
  python excess.py --run                          # all defaults: 50% excess,
                                                  #   7-yr linear grade-out, 2025
  python excess.py --excess-all 60 --grade-out 10
  python excess.py --valuation-year 2027 --age 65 --sex male

----------------------------------------------------------------------------
DRIVING IT FROM THE PULLFORWARD SIDE INSTEAD
----------------------------------------------------------------------------
  python excess.py --peak 0.65 --grade-out 7      # assume 65% of 2021's deaths
                                                  #   were pulled into 2020; the
                                                  #   implied excess % is reported
  python excess.py --peak 0.5 --pullforward-shape exponential --decay-rate 0.3

----------------------------------------------------------------------------
HOW THE PULLFORWARD (HARVEST) GRADES AWAY
----------------------------------------------------------------------------
  python excess.py --pullforward-shape linear        # ramps to zero at --grade-out (default)
  python excess.py --pullforward-shape step           # full effect inside --grade-out, then zero
  python excess.py --pullforward-shape exponential --decay-rate 0.3   # fat tail, no cutoff

----------------------------------------------------------------------------
HOW THE EXCESS ITSELF IS TIMED ACROSS YEARS
----------------------------------------------------------------------------
  python excess.py --excess-shape instant              # all in 2020 (default)
  python excess.py --excess-shape linear --excess-spread 3   # fades to zero over 3 yrs
  python excess.py --excess-shape exponential --excess-decay 0.5
                                                        #   x(j) ~ e^(-r*j): levels
                                                        #   100/61/37/22/...% of 2020's,
                                                        #   halving every ln2/r yrs; cut
                                                        #   below 0.1% and renormalized
  python excess.py --excess-shape gaussian              # Gaussian fade w(j) = 2^(-j^2/4),
                                                        #   fitted to the actual COVID
                                                        #   pattern: 100/84/50/21/6% of the
                                                        #   2020 level over 2020-2024, zero
                                                        #   from 2025 (38/32/19/8/2% of the
                                                        #   total excess). "empirical" is
                                                        #   accepted as an alias

----------------------------------------------------------------------------
PER-BAND EXCESS (21 five-year bands: 0-4, 5-9, ..., 95-99, 100+)
----------------------------------------------------------------------------
  python excess.py --excess "30,30,30,30,30,35,35,40,40,45,45,50,50,55,55,60,60,60,55,50,45"

  # sex-specific bands (either flag overrides --excess/--excess-all for that sex):
  python excess.py --excess-male "60,...,55" --excess-female "45,...,40"

----------------------------------------------------------------------------
MORTALITY IMPROVEMENT (same flags as sensitivity.py; default flat 1%/yr)
----------------------------------------------------------------------------
  python excess.py --run --improvement-rate 0.015
  python excess.py --run --no-improvement
  python excess.py --run --improvement-table my_scale.csv

----------------------------------------------------------------------------
OUTPUT
----------------------------------------------------------------------------
  --age / --sex        focus on one cohort (ages are the age IN 2020)
  --trajectory-age N   age for the printed 2010-2035 trajectory (default 65)
  --save               write CSVs to output/excess_calibration/

Notes
-----
* Harvested deaths come from the cohort AS IT AGES: a 60-year-old's death
  pulled forward from 2023 would have happened at age 63. Mortality rises
  with age, so a fixed number of harvested deaths is a shrinking share of
  each later year's deaths — the solved pullforward percentages are much
  smaller than the headline excess percentage. This holds for every
  --pullforward-shape: `peak` is always "share of 2021's deaths pulled into
  2020", since every shape is normalized to 1.0 at t=1.
* --excess-shape only re-times WHEN the excess deaths happened (the totals
  and the solved pullforward peak are unchanged); --pullforward-shape governs
  how the harvesting (the repayment) grades away over time.
"""
from __future__ import annotations

import argparse
import os
import sys

import pandas as pd

from mortality_model.ssa_table import load_ssa_table
from mortality_model.improvement import load_improvement_scale
from mortality_model.excess import (
    N_BANDS,
    band_label,
    default_excess_bands,
    run_excess_cohort,
    mortality_trajectory,
    TRAJECTORY_START,
    TRAJECTORY_END,
    PULLFORWARD_SHAPES,
    EXCESS_SHAPES,
    DEFAULT_DECAY_RATE,
    DEFAULT_EXCESS_DECAY,
    exp_excess_years,
)

DEFAULT_AGES = [40, 50, 55, 60, 65, 70, 75, 80]  # ages IN 2020


def _parse_band_list(text: str, flag: str) -> list[float]:
    parts = [p for p in text.replace(" ", "").split(",") if p != ""]
    if len(parts) != N_BANDS:
        bands = ", ".join(band_label(i) for i in range(N_BANDS))
        print(f"{flag} needs exactly {N_BANDS} comma-separated percentages, "
              f"one per band:\n  {bands}\nGot {len(parts)} values.")
        sys.exit(1)
    vals = [float(p) for p in parts]
    for v in vals:
        if not (0.0 <= v <= 1000.0):
            print(f"Excess percentages must be between 0 and 1000 (got {v}).")
            sys.exit(1)
    return vals


def parse_excess(args) -> dict:
    """Per-sex 21-band excess vectors (fractions): {'male': [...], 'female': [...]}."""
    if args.excess is not None:
        base = _parse_band_list(args.excess, "--excess")
    else:
        if not (0.0 <= args.excess_all <= 1000.0):
            print(f"--excess-all must be between 0 and 1000 (got {args.excess_all}).")
            sys.exit(1)
        base = [args.excess_all] * N_BANDS
    out = {"male": base, "female": list(base)}
    if args.excess_male is not None:
        out["male"] = _parse_band_list(args.excess_male, "--excess-male")
    if args.excess_female is not None:
        out["female"] = _parse_band_list(args.excess_female, "--excess-female")
    return {sx: [v / 100.0 for v in vals] for sx, vals in out.items()}


def build_scale(args):
    if args.no_improvement:
        return None
    rate = args.improvement_rate if args.improvement_rate is not None else 0.01
    return load_improvement_scale(args.improvement_table, rate)


def describe(args, excess_bands) -> str:
    if args.pullforward_shape == "exponential":
        harvest_desc = f"exponential decay, rate {args.decay_rate} (no hard cutoff)"
    elif args.pullforward_shape == "step":
        harvest_desc = f"step: full effect for {args.grade_out} years, then zero"
    else:
        harvest_desc = f"linear: ramps to zero over {args.grade_out} years"
    if args.excess_shape == "linear":
        excess_desc = f"graded linearly to zero over {args.excess_spread} years"
    elif args.excess_shape == "exponential":
        yrs = exp_excess_years(args.excess_decay)
        excess_desc = (f"exponential fade x(j) ~ e^(-{args.excess_decay:g}*j): "
                       f"halves every {0.6931 / args.excess_decay:.1f} yrs, "
                       f"cut below 0.1% of 2020's level ({yrs} years)")
    elif args.excess_shape == "gaussian":
        excess_desc = ("Gaussian fade w(j) = 2^(-j^2/4), fitted to the observed pattern: "
                       "100/84/50/21/6% of the 2020 level over 2020-2024, zero after")
    else:
        excess_desc = "all in 2020"
    if args.peak is not None:
        driver_line = (f"  Pullforward peak  : {args.peak * 100:.0f}% of 2021's deaths "
                       "(ASSUMED; the cumulative excess is implied)")
    else:
        same = excess_bands["male"] == excess_bands["female"]
        uniform = same and len(set(excess_bands["male"])) == 1
        driver_line = ("  Cumulative excess : "
                       + (f"{excess_bands['male'][0] * 100:.0f}% of one year's deaths, all ages/sexes"
                          if uniform
                          else ("per-band (see --excess)" if same
                                else "per-band and per-sex (see --excess-male / --excess-female)"))
                       + " (ASSUMED; the pullforward peak is solved)")
    lines = [
        driver_line,
        f"  Harvest shape     : {harvest_desc}",
        f"  Excess timing     : {excess_desc}",
        f"  Valuation year    : {args.valuation_year}",
        f"  Improvement       : "
        + ("OFF (static 2019 table)" if args.no_improvement
           else (f"from CSV: {args.improvement_table}" if args.improvement_table
                 else f"flat {(args.improvement_rate if args.improvement_rate is not None else 0.01):.2%}/yr, "
                      f"anchored to the 2019 table")),
    ]
    return "\n".join(lines)


def main():
    p = argparse.ArgumentParser(
        description="COVID excess-mortality calibration mode.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--run", action="store_true",
                   help="run with all defaults (used when no other flag is given)")
    p.add_argument("--excess-all", type=float, default=50.0, dest="excess_all",
                   help="cumulative excess for EVERY band, %% of one-year deaths (default 50)")
    p.add_argument("--excess",
                   help=f"{N_BANDS} comma-separated percentages, one per 5-year band "
                        "(0-4, 5-9, ..., 95-99, 100+)")
    p.add_argument("--excess-male", dest="excess_male",
                   help=f"{N_BANDS} comma-separated percentages for MALES only "
                        "(overrides --excess/--excess-all for that sex)")
    p.add_argument("--excess-female", dest="excess_female",
                   help=f"{N_BANDS} comma-separated percentages for FEMALES only "
                        "(overrides --excess/--excess-all for that sex)")
    p.add_argument("--peak", type=float,
                   help="drive the model from the pullforward side instead: the ASSUMED "
                        "share of 2021's deaths pulled into 2020 (0..1). Overrides "
                        "--excess/--excess-all; the implied cumulative excess is reported")
    p.add_argument("--grade-out", type=int, default=7, dest="grade_out",
                   help="years over which the harvesting grades away, for --pullforward-shape "
                        "linear/step (default 7; unused for exponential)")
    p.add_argument("--pullforward-shape", choices=PULLFORWARD_SHAPES, default="linear",
                   dest="pullforward_shape",
                   help="how the harvest (repayment) grades away: linear (default), "
                        "step, or exponential")
    p.add_argument("--decay-rate", type=float, default=DEFAULT_DECAY_RATE, dest="decay_rate",
                   help=f"exponential decay rate, only used with --pullforward-shape "
                        f"exponential (default {DEFAULT_DECAY_RATE})")
    p.add_argument("--valuation-year", type=int, default=2025, dest="valuation_year",
                   help="standing point for LE gain / remaining multiples (default 2025)")
    p.add_argument("--excess-shape", choices=EXCESS_SHAPES + ("empirical",),
                   default="instant", dest="excess_shape",
                   help="how the 2020 excess itself is timed: instant (default, all in "
                        "2020), linear (fades over --excess-spread years), exponential "
                        "(x(j) ~ e^(-r*j) with r = --excess-decay), or gaussian "
                        "(fade 2^(-j^2/4) fitted to the observed COVID pattern: "
                        "100/84/50/21/6%% of the 2020 level over 2020-2024, zero from "
                        "2025; 'empirical' is accepted as an alias)")
    p.add_argument("--excess-spread", type=int, default=3, dest="excess_spread",
                   help="years over which the excess fades to zero, for --excess-shape "
                        "linear (default 3)")
    p.add_argument("--excess-decay", type=float, default=DEFAULT_EXCESS_DECAY,
                   dest="excess_decay",
                   help=f"decay rate for --excess-shape exponential (default "
                        f"{DEFAULT_EXCESS_DECAY}: the excess halves every "
                        f"~{0.6931 / DEFAULT_EXCESS_DECAY:.1f} years)")
    p.add_argument("--age", type=int, help="focus on one age IN 2020")
    p.add_argument("--sex", choices=["male", "female"])
    p.add_argument("--trajectory-age", type=int, dest="trajectory_age",
                   help="age for the 2010-2035 trajectory table (default: --age or 65)")
    p.add_argument("--improvement-rate", type=float, dest="improvement_rate",
                   help="flat annual improvement rate (0.01 = 1%%/yr; default 0.01)")
    p.add_argument("--improvement-table", dest="improvement_table",
                   help="path to a 1D or 2D improvement CSV")
    p.add_argument("--no-improvement", action="store_true", dest="no_improvement")
    p.add_argument("--save", action="store_true",
                   help="write CSVs to output/excess_calibration/")

    if len(sys.argv) == 1:
        print(__doc__)
        return
    args = p.parse_args()

    if args.peak is not None and not (0.0 <= args.peak <= 1.0):
        print("--peak must be between 0 and 1 (e.g. 0.65 = 65% of 2021's deaths).")
        sys.exit(1)
    if args.pullforward_shape != "exponential" and args.grade_out < 1:
        print("--grade-out must be at least 1 year (for linear/step --pullforward-shape).")
        sys.exit(1)
    if args.pullforward_shape == "exponential" and args.decay_rate <= 0:
        print("--decay-rate must be positive.")
        sys.exit(1)
    if args.excess_shape == "empirical":
        args.excess_shape = "gaussian"
    if args.excess_shape == "linear" and args.excess_spread < 2:
        print("--excess-spread must be at least 2 years (for --excess-shape linear).")
        sys.exit(1)
    if args.excess_shape == "exponential" and args.excess_decay <= 0:
        print("--excess-decay must be positive.")
        sys.exit(1)

    excess_bands = parse_excess(args)
    table = load_ssa_table(None, 119)
    scale = build_scale(args)

    ages = [args.age] if args.age is not None else DEFAULT_AGES
    sexes = [args.sex] if args.sex is not None else ["male", "female"]

    print("=" * 78)
    print("COVID pullforward model (calendar-anchored)")
    print(describe(args, excess_bands))
    print("=" * 78)

    results = []
    for sex in sexes:
        for age in ages:
            results.append(run_excess_cohort(
                age, sex, table, scale, excess_bands[sex], args.grade_out,
                args.valuation_year, args.pullforward_shape, args.decay_rate,
                args.excess_shape, args.excess_spread, 119, args.peak,
                args.excess_decay,
            ))

    direct = args.peak is not None
    excess_col = "Implied excess (%)" if direct else "Excess (%)"
    peak_col = "Pulled from 2021 (%)" + (" [given]" if direct else " [solved]")
    rows = []
    any_infeasible = False
    for r in results:
        any_infeasible = any_infeasible or r["infeasible"]
        rows.append({
            "Sex": r["sex"].capitalize(),
            "Age 2020": r["age_2020"],
            f"Age {args.valuation_year}": r["age_at_valuation"],
            excess_col: round(r["excess_fraction"] * 100, 1),
            peak_col: round(r["peak"] * 100, 2),
            "Alive vs baseline (%)": round(r["alive_vs_baseline"] * 100, 2),
            "LE base (yrs)": round(r["le_base"], 2),
            "LE surv (yrs)": round(r["le_surv"], 2),
            "LE gain (yrs)": round(r["le_change"], 3),
            "Equiv mult (%)": round(r["equiv_mult"] * 100, 2),
            "!": "INFEASIBLE" if r["infeasible"] else "",
        })
    df = pd.DataFrame(rows).set_index(["Sex", "Age 2020"])
    print(f"\n=== Standing in {args.valuation_year}: survivors vs the no-COVID baseline ===")
    print(df.to_string())
    if direct:
        print("\n  Implied excess (%): the cumulative excess deaths (as % of one normal")
        print("  year's deaths) that the ASSUMED pullforward implies, via the same")
        print("  conservation equation the excess driver solves in reverse.")
    else:
        print("\n  Pulled from 2021 (%): the SOLVED peak pullforward -- the share of 2021's")
        print("  deaths that instead happened in 2020. It is much smaller than the excess %")
        print("  because the harvested deaths spread over the horizon AND the cohort ages")
        print("  into higher mortality, so each borrowed death is a smaller share of that")
        print("  year's (larger) death count.")
    print("  Alive vs baseline (%): population still alive at the valuation year,")
    print("  relative to the no-COVID baseline (the not-yet-harvested deficit).")
    print("  Equiv mult (%) -- THE KEY OUTPUT: the single flat qx multiplier on the")
    print("  whole baseline table (from the valuation year on) that reproduces the")
    print("  survivors' LE, e.g. 99.5 = a uniform 0.5% mortality cut at every age.")
    if any_infeasible:
        fix = ("increase --decay-rate's inverse (use a smaller --decay-rate for a fatter tail)"
               if args.pullforward_shape == "exponential" else "lengthen --grade-out")
        print("\n  WARNING - INFEASIBLE rows: the assumed excess exceeds ALL deaths available")
        print("  to harvest under this shape (the solved pullforward would exceed 100% of")
        print(f"  2021's deaths). The pullforward was capped at 100%; {fix}")
        print("  or lower the excess for those ages.")

    # Trajectory table
    traj_age = args.trajectory_age if args.trajectory_age is not None \
        else (args.age if args.age is not None else 65)
    trajectories = []
    for sex in sexes:
        traj = mortality_trajectory(
            traj_age, sex, table, scale, excess_bands[sex], args.grade_out,
            args.valuation_year, args.pullforward_shape, args.decay_rate,
            args.excess_shape, args.excess_spread, peak=args.peak,
            excess_decay_rate=args.excess_decay,
        )
        trajectories.append(traj)
        print(f"\n=== Mortality trajectory at age {traj_age}, {sex} "
              f"({TRAJECTORY_START}-{TRAJECTORY_END}) ===")
        print(f"  {'Year':>6}  {'Baseline qx':>12}  {'With COVID':>12}  {'Ratio':>7}")
        for i, y in enumerate(traj["years"]):
            b, c = traj["baseline_qx"][i], traj["covid_qx"][i]
            ratio = c / b if b > 0 else float("nan")
            marker = "  <- excess" if ratio > 1.001 else (" <- harvest" if ratio < 0.999 else "")
            print(f"  {y:>6}  {b:>12.5f}  {c:>12.5f}  {ratio:>7.3f}{marker}")

    if args.save:
        folder = os.path.join("output", "excess_calibration")
        os.makedirs(folder, exist_ok=True)
        summary = pd.DataFrame([{
            "sex": r["sex"],
            "age_2020": r["age_2020"],
            "age_at_valuation": r["age_at_valuation"],
            "driver": r["driver"],
            "excess_pct": r["excess_fraction"] * 100,
            "pull_from_2021_pct": r["peak"] * 100,
            "infeasible": r["infeasible"],
            "alive_vs_baseline": r["alive_vs_baseline"],
            "le_base": r["le_base"],
            "le_surv": r["le_surv"],
            "le_gain": r["le_change"],
            "equiv_flat_qx_multiplier": r["equiv_mult"],
        } for r in results])
        path = os.path.join(folder, "summary.csv")
        summary.to_csv(path, index=False)
        print(f"\n  Saved: {path}")
        for r in results:
            detail = pd.DataFrame({
                "year": r["years"], "age": r["ages"],
                "baseline_qx": r["q_base"], "covid_qx": r["q_covid"],
                "multiple": r["multiple"], "pullforward_f": r["f"],
                "excess_deaths": r["x_excess"], "harvested_deaths": r["h_harvest"],
                "baseline_deaths": r["D_base"], "covid_deaths": r["D_covid"],
                "baseline_alive": r["A_base"][:-1], "covid_alive": r["A_covid"][:-1],
            })
            detail.to_csv(os.path.join(
                folder, f"path_{r['sex']}_age{r['age_2020']}.csv"), index=False)
        for traj in trajectories:
            pd.DataFrame({
                "year": traj["years"],
                "baseline_qx": traj["baseline_qx"],
                "covid_qx": traj["covid_qx"],
            }).to_csv(os.path.join(
                folder, f"trajectory_{traj['sex']}_age{traj['age']}.csv"), index=False)
        print(f"  Saved: {folder}/path_*.csv, trajectory_*.csv")


if __name__ == "__main__":
    main()
