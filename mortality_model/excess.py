"""
COVID excess-mortality calibration mode.

Instead of assuming the pullforward curve f(t) directly, this mode starts from
an observable: the CUMULATIVE excess mortality each age group suffered because
of COVID, expressed as a percentage of that group's one-year death rate in
2020. Example: excess of 50% for 60-year-olds means the group's extra COVID
deaths equalled half of a normal year's deaths for 60-year-olds.

All excess deaths are attributed to a single year (2020). Because the people
who died early were "borrowed" from future years, the same number of deaths
must be MISSING from 2021 onward — and, crucially, those missing deaths come
from the same cohort as it ages (a 60-year-old's pulled-forward 2023 death
would have happened at age 63, not 60). Since mortality rises with age, a
fixed number of harvested deaths is a shrinking share of each later year's
deaths — the aging of the cohort automatically mutes the future effect.

Given the excess E and a pullforward (harvest) shape, the model SOLVES for
the one remaining unknown: the peak share of next year's (2021's) deaths that
were pulled into 2020.

    excess deaths            harvested deaths
    E x qx(age, 2020)   =    sum over t>=1 of  peak x w(t) x D_b(t)

where w(t) is the chosen HARVEST SHAPE (normalized so w(1) = 1.0, i.e. `peak`
is always exactly "the share of 2021 deaths pulled into 2020" regardless of
shape) and D_b(t) is the cohort's baseline (no-COVID) unconditional deaths in
year 2020+t:

    linear:      w(t) = max(0, 1 - (t-1)/G)              for t = 1..G, else 0
    step:        w(t) = 1                                for t = 1..G, else 0
    exponential: w(t) = exp(-decay_rate x (t-1))          for all t >= 1 (no cutoff)

Everything is bookkept in absolute deaths per person alive at the start of
2020:

    baseline:  D_b(t) = A_b(t) x qx(age+t, 2020+t),  A_b(t+1) = A_b(t) - D_b(t)
    COVID:     D_c(t) = D_b(t) + x(t) - h(t),        A_c(t+1) = A_c(t) - D_c(t)

with x(t) the TIMING of the excess deaths themselves and h(t) = f(t) x D_b(t)
the harvested deaths. x(t) has its own shape, independent of the harvest
shape:

    instant:   all of X lands in 2020 (t=0)
    linear:    X fades linearly to zero over N years from 2020 (a triangular
               spread: 2020 gets the most, year N-1 the least)
    empirical: a preset shape matching the actual 2020-2022 COVID experience —
               2020 and 2021 near their peak (2021 slightly lower), about half
               that level in 2022, a muted tail in 2023, almost nothing left
               by 2024, and zero from 2025 on (EMPIRICAL_EXCESS_WEIGHTS)

Period death rates, mortality multiples, and the valuation-year life
expectancies follow directly:

    q_c(t) = D_c(t) / A_c(t),   multiple(t) = q_c(t) / q_b(t)
    LE at valuation year V = sum over k>=0 of A(tv+k) / A(tv),  tv = V - 2020

For linear/step harvest shapes the multiple returns to exactly 1.0 once the
grade-out horizon (and any later-ending excess timing window) has passed; for
exponential harvest there is no hard horizon, so it only approaches 1.0.

Calendar anchoring: the base table is SSA 2019, so qx(age, year) applies the
mortality improvement scale forward from 2019 (and BACKWARD for years before
2019, by dividing out the improvement) — this differs from the direct
pullforward mode, which counts projection years from 2022.
"""
from __future__ import annotations

import math
from typing import Optional

import numpy as np
import pandas as pd

from mortality_model.ssa_table import get_qx
from mortality_model.improvement import ImprovementScale
from mortality_model.analysis import _le_from_qx

TABLE_BASE_YEAR = 2019   # calendar year of the SSA 2019 base table
COVID_YEAR = 2020        # the year all excess deaths are attributed to
N_BANDS = 21             # 5-year bands: 0-4, 5-9, ..., 95-99, 100+
DEFAULT_EXCESS = 0.50    # cumulative excess as a fraction of one-year mortality
TRAJECTORY_START = 2010
TRAJECTORY_END = 2035
DEFAULT_DECAY_RATE = 0.3

# Relative excess "level" by year, 2020..2024 (index 0..4); zero from 2025 on.
# Matches the empirical COVID pattern: 2020 and 2021 near their peak (2021
# slightly lower), about half that level in 2022, a muted tail in 2023,
# almost nothing left by 2024. Normalized to sum to 1 before scaling by X.
EMPIRICAL_EXCESS_WEIGHTS = (1.00, 0.90, 0.50, 0.20, 0.05)

PULLFORWARD_SHAPES = ("linear", "step", "exponential")
EXCESS_SHAPES = ("instant", "linear", "empirical")


def band_index(age: int) -> int:
    """Index of the 5-year excess band containing `age` (last band = 100+)."""
    return min(max(age, 0) // 5, N_BANDS - 1)


def band_label(i: int) -> str:
    return "100+" if i == N_BANDS - 1 else f"{5 * i}-{5 * i + 4}"


def default_excess_bands(excess: float = DEFAULT_EXCESS) -> list[float]:
    """One cumulative-excess fraction per 5-year band (decimals: 0.5 = 50%)."""
    return [excess] * N_BANDS


def calendar_factor(scale: Optional[ImprovementScale], age: int, year: int,
                    sex: str) -> float:
    """
    Cumulative improvement factor taking the base table (2019) to `year`.
    Years after 2019 compound (1 - rate); years before 2019 divide it out,
    so the same scale also back-projects the table to 2010.
    """
    if scale is None or year == TABLE_BASE_YEAR:
        return 1.0
    if scale.mode in ("flat", "1d"):
        # Rate does not depend on the calendar year; a negative exponent
        # de-improves the table for years before 2019.
        r = scale.rate(age, year, sex)
        return (1.0 - r) ** (year - TABLE_BASE_YEAR)
    if year > TABLE_BASE_YEAR:
        f = 1.0
        for k in range(TABLE_BASE_YEAR + 1, year + 1):
            f *= 1.0 - scale.rate(age, k, sex)
        return f
    f = 1.0
    for k in range(year + 1, TABLE_BASE_YEAR + 1):
        f *= 1.0 - scale.rate(age, k, sex)
    return 1.0 / f


def calendar_qx(table: pd.DataFrame, age: int, year: int, sex: str,
                scale: Optional[ImprovementScale], max_age: int = 119) -> float:
    """qx for (age, calendar year): the 2019 table +/- improvement."""
    a = min(age, max_age)
    if a >= max_age:
        return 1.0
    return min(1.0, get_qx(table, a, sex) * calendar_factor(scale, a, year, sex))


def harvest_weight(t: int, shape: str, grade_out_years: int,
                   decay_rate: float = DEFAULT_DECAY_RATE) -> float:
    """
    w(t): un-scaled harvest weight for projection year t (1-indexed calendar
    offset from 2020; t <= 0 is never harvested). Always normalized so
    w(1) == 1.0, so the solved `peak` is always exactly "share of 2021's
    deaths pulled into 2020" regardless of shape.
    """
    if t <= 0:
        return 0.0
    if shape == "step":
        return 1.0 if (grade_out_years > 0 and t <= grade_out_years) else 0.0
    if shape == "exponential":
        return math.exp(-decay_rate * (t - 1))
    # linear (default)
    if grade_out_years <= 0:
        return 0.0
    return max(0.0, 1.0 - (t - 1) / grade_out_years)


def excess_timing(X: float, shape: str, spread_years: int, n: int) -> list[float]:
    """
    x(t): absolute excess deaths landing in projection year t (0-indexed,
    t=0 is 2020), for the chosen timing shape. Always sums to X (up to
    truncation if the cohort's remaining lifetime n is very short).
    """
    x = [0.0] * n
    if shape == "linear" and spread_years and spread_years >= 2:
        gx = int(spread_years)
        wsum = gx * (gx + 1) / 2
        for j in range(min(gx, n)):
            x[j] = X * (gx - j) / wsum
    elif shape == "empirical":
        w = EMPIRICAL_EXCESS_WEIGHTS
        wsum = sum(w)
        for j in range(min(len(w), n)):
            x[j] = X * w[j] / wsum
    else:  # "instant"
        x[0] = X
    return x


def excess_window_end(shape: str, spread_years: int) -> int:
    """Last t-index (0-indexed) with potentially nonzero excess timing."""
    if shape == "linear" and spread_years and spread_years >= 2:
        return int(spread_years) - 1
    if shape == "empirical":
        return len(EMPIRICAL_EXCESS_WEIGHTS) - 1
    return 0


def _solve_equiv_multiplier(target_le: float, qx_vec: np.ndarray) -> float:
    """Flat qx multiplier on the future baseline that reproduces target_le."""
    lo, hi = 1e-6, 1.0
    while _le_from_qx(qx_vec * hi) > target_le and hi < 64.0:
        lo = hi
        hi *= 2.0
    if _le_from_qx(qx_vec * lo) < target_le:
        return float("nan")
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if _le_from_qx(qx_vec * mid) > target_le:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def run_excess_cohort(
    age_2020: int,
    sex: str,
    table: pd.DataFrame,
    scale: Optional[ImprovementScale],
    excess_by_band: list[float],
    grade_out_years: int,
    valuation_year: int = 2025,
    pullforward_shape: str = "linear",
    decay_rate: float = DEFAULT_DECAY_RATE,
    excess_shape: str = "instant",
    excess_spread_years: int = 3,
    max_age: int = 119,
) -> dict:
    """
    Full excess-calibration bookkeeping for one cohort (aged `age_2020` in
    2020). Index t = 0, 1, 2, ... is calendar year 2020 + t; the cohort is
    aged age_2020 + t in that year. All quantities are per person alive at
    the start of 2020.
    """
    n = max_age - age_2020 + 1

    # Baseline path: calendar-anchored qx along the cohort's diagonal
    q = [calendar_qx(table, age_2020 + t, COVID_YEAR + t, sex, scale, max_age)
         for t in range(n)]
    q[n - 1] = 1.0
    A_b = [0.0] * (n + 1)
    D_b = [0.0] * n
    A_b[0] = 1.0
    for t in range(n):
        D_b[t] = A_b[t] * q[t]
        A_b[t + 1] = A_b[t] - D_b[t]

    # Excess deaths in 2020: E x one year's mortality
    E = excess_by_band[band_index(age_2020)]
    X = E * q[0]

    # Solve the pullforward peak so harvested deaths (2021+) equal the excess,
    # under the chosen harvest shape (w(1) is always normalized to 1.0).
    G = grade_out_years
    w = [harvest_weight(t, pullforward_shape, G, decay_rate) for t in range(n)]
    denom = 0.0
    for t in range(n):
        denom += w[t] * D_b[t]
    peak = X / denom if denom > 0 else math.inf
    infeasible = (not math.isfinite(peak)) or peak > 1.0 + 1e-12

    f = [min(1.0, peak * w[t]) if w[t] > 0 else 0.0 for t in range(n)]
    h = [f[t] * D_b[t] for t in range(n)]
    harvested = 0.0
    for t in range(n):
        harvested += h[t]
    shortfall = X - harvested  # > 0 only when infeasible (f capped at 100%)

    # Timing of the excess deaths, under the chosen excess-timing shape.
    x = excess_timing(X, excess_shape, excess_spread_years, n)

    # COVID path. Inside the excess/harvest window, deaths are baseline deaths
    # plus the excess timing minus the harvest (absolute-death bookkeeping).
    # After the window, survivors revert to baseline RATES: identical to the
    # baseline path whenever the harvest fully repaid the excess, and the
    # graceful behavior when it could not (infeasible inputs). Exponential
    # harvest has no hard horizon, so it never reverts -- the raw formula is
    # used throughout, and h(t)/x(t) simply decay toward (numerically
    # indistinguishable from) zero.
    if pullforward_shape in ("linear", "step"):
        window = max(G, excess_window_end(excess_shape, excess_spread_years))
    else:
        window = n  # sentinel: "elif t > window" below never fires
    A_c = [0.0] * (n + 1)
    D_c = [0.0] * n
    A_c[0] = 1.0
    for t in range(n):
        if t == n - 1:
            dc = A_c[t]  # terminal age: qx = 1 regardless of harvest
        elif t > window:
            dc = A_c[t] * q[t]
        else:
            dc = D_b[t] + x[t] - h[t]
        dc = max(0.0, min(dc, A_c[t]))
        D_c[t] = dc
        A_c[t + 1] = A_c[t] - dc

    q_c = [D_c[t] / A_c[t] if A_c[t] > 1e-15 else float("nan") for t in range(n)]
    multiple = [q_c[t] / q[t] if q[t] > 0 and math.isfinite(q_c[t]) else float("nan")
                for t in range(n)]

    # Valuation: condition on being alive at the start of valuation_year
    tv = valuation_year - COVID_YEAR
    valid = 0 <= tv < n and A_c[tv] > 1e-12 and A_b[tv] > 1e-12
    le_base = le_surv = le_change = equiv_mult = alive_vs_baseline = float("nan")
    if valid:
        le_base = 0.0
        for k in range(tv, n):
            le_base += A_b[k] / A_b[tv]
        le_surv = 0.0
        for k in range(tv, n):
            le_surv += A_c[k] / A_c[tv]
        le_change = le_surv - le_base
        equiv_mult = _solve_equiv_multiplier(le_surv, np.array(q[tv:]))
        alive_vs_baseline = A_c[tv] / A_b[tv]

    return {
        "age_2020": age_2020,
        "sex": sex,
        "n": n,
        "years": [COVID_YEAR + t for t in range(n)],
        "ages": [age_2020 + t for t in range(n)],
        "q_base": q,
        "D_base": D_b,
        "A_base": A_b,
        "q_covid": q_c,
        "D_covid": D_c,
        "A_covid": A_c,
        "f": f,
        "x_excess": x,
        "h_harvest": h,
        "multiple": multiple,
        "excess_fraction": E,
        "excess_deaths": X,
        "peak": peak,                    # solved f(2021): share of 2021 deaths pulled into 2020
        "infeasible": infeasible,
        "harvested": harvested,
        "shortfall": shortfall,
        "valuation_year": valuation_year,
        "tv": tv,
        "valid": valid,
        "age_at_valuation": age_2020 + tv,
        "le_base": le_base,
        "le_surv": le_surv,
        "le_change": le_change,
        "equiv_mult": equiv_mult,
        "alive_vs_baseline": alive_vs_baseline,
    }


def mortality_trajectory(
    age: int,
    sex: str,
    table: pd.DataFrame,
    scale: Optional[ImprovementScale],
    excess_by_band: list[float],
    grade_out_years: int,
    valuation_year: int = 2025,
    pullforward_shape: str = "linear",
    decay_rate: float = DEFAULT_DECAY_RATE,
    excess_shape: str = "instant",
    excess_spread_years: int = 3,
    year_start: int = TRAJECTORY_START,
    year_end: int = TRAJECTORY_END,
    max_age: int = 119,
    cohort_cache: Optional[dict] = None,
) -> dict:
    """
    Period mortality for a FIXED age across calendar years: the baseline
    trajectory (table +/- improvement) and the COVID-impacted trajectory.
    The COVID rate for (age, 2020+k) comes from the cohort that was age-k in
    2020 — each year's rate reflects a different cohort's spike/harvest.
    """
    years = list(range(year_start, year_end + 1))
    base, covid = [], []
    cache = {} if cohort_cache is None else cohort_cache
    for y in years:
        b = calendar_qx(table, age, y, sex, scale, max_age)
        base.append(b)
        k = y - COVID_YEAR
        c = age - k
        if k < 0 or c < 0 or age > max_age:
            covid.append(b)
            continue
        key = (c, sex)
        if key not in cache:
            cache[key] = run_excess_cohort(
                c, sex, table, scale, excess_by_band, grade_out_years,
                valuation_year, pullforward_shape, decay_rate,
                excess_shape, excess_spread_years, max_age,
            )
        r = cache[key]
        covid.append(r["q_covid"][k] if k < r["n"] else float("nan"))
    return {"age": age, "sex": sex, "years": years,
            "baseline_qx": base, "covid_qx": covid}
