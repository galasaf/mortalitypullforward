"""
Mortality improvement scales.

A mortality improvement scale gives MI(age, year, sex): the annual rate at
which qx improves (declines). Improvement compounds from the base year
(cfg.covid_end_year) forward:

    qx(age, year) = qx_base(age) * PRODUCT over k in (base_year, year] of (1 - MI(age, k, sex))

Projection year t (1-indexed) corresponds to calendar year base_year + t - 1,
so year 1 of the projection uses the base table unchanged, year 2 gets one
year of improvement, and so on. With a flat 1% scale the year-t qx is
qx_base * 0.99^(t-1).

Three kinds of scale are supported:
  flat — one rate for every age, year, and sex (the default: 1%)
  1D   — a rate per age (constant across calendar years), optionally by sex
  2D   — a rate per age AND calendar year (MP-style), optionally by sex

Scales are imported from CSV (Excel: "Save As -> CSV"). The format is
auto-detected from the header:

1D CSV schema (no year columns):
    age,improvement                 <- one rate applied to both sexes, OR
    age,male_improvement,female_improvement
    0,0.010
    1,0.010
    ...

2D CSV schema (wide; every column after age/sex is a 4-digit calendar year):
    age,2023,2024,2025,...          <- applies to both sexes, OR
    age,sex,2023,2024,2025,...      <- sex is 'male' or 'female'
    0,0.010,0.010,0.010,...

Rules for both schemas:
  * Rates are DECIMALS: 0.01 means 1% improvement per year. Negative values
    (mortality deterioration) are allowed. Values with |rate| > 0.20 are
    rejected — they almost always mean the file used percents (1 vs 0.01).
  * Ages outside the range in the file use the nearest edge age.
  * (2D) Calendar years outside the range in the file use the nearest edge
    year — e.g. years after the last column hold that column's rate forever.
  * (2D) If a 'sex' column is present, each sex uses its own rows; a sex with
    no rows gets zero improvement.

Use write_template_1d() / write_template_2d() (or
`python sensitivity.py --improvement-template 1d|2d`) to generate a
pre-filled template to edit.
"""
from __future__ import annotations

import re
from typing import Optional

import numpy as np
import pandas as pd

_YEAR_RE = re.compile(r"^\d{4}$")
_MAX_SANE_RATE = 0.20


class ImprovementScale:
    """
    Callable mortality improvement scale: rate(age, year, sex) -> float.

    Internally one of:
      mode='flat' : self.flat_rate
      mode='1d'   : self.rates_1d[sex] = np.ndarray indexed by age 0..max_age
      mode='2d'   : self.rates_2d[sex] = np.ndarray [max_age+1, n_years],
                    self.years = sorted list of calendar years (columns)
    """

    def __init__(self, mode: str, flat_rate: float = 0.0,
                 rates_1d: dict = None, rates_2d: dict = None,
                 years: list = None, source: str = None):
        self.mode = mode
        self.flat_rate = flat_rate
        self.rates_1d = rates_1d or {}
        self.rates_2d = rates_2d or {}
        self.years = years or []
        self.source = source  # file path, or None for flat

    def rate(self, age: int, year: int, sex: str) -> float:
        """Annual improvement rate for (age, calendar year, sex)."""
        if self.mode == "flat":
            return self.flat_rate

        if self.mode == "1d":
            arr = self.rates_1d.get(sex)
            if arr is None:
                return 0.0
            return float(arr[min(max(age, 0), len(arr) - 1)])

        # 2d
        arr = self.rates_2d.get(sex)
        if arr is None:
            return 0.0
        a = min(max(age, 0), arr.shape[0] - 1)
        y = min(max(year, self.years[0]), self.years[-1])
        return float(arr[a, self.years.index(y)])

    def cumulative_factor(self, age: int, base_year: int, n_years: int, sex: str) -> float:
        """
        PRODUCT over j=1..n_years of (1 - rate(age, base_year + j, sex)).

        This is the total qx adjustment for someone reaching `age` after
        `n_years` of improvement from the base year. n_years=0 returns 1.0.
        """
        if n_years <= 0:
            return 1.0
        if self.mode in ("flat", "1d"):
            # Rate does not depend on the calendar year.
            r = self.rate(age, base_year + 1, sex)
            return (1.0 - r) ** n_years
        factor = 1.0
        for j in range(1, n_years + 1):
            factor *= 1.0 - self.rate(age, base_year + j, sex)
        return factor

    def describe(self) -> str:
        if self.mode == "flat":
            return f"flat {self.flat_rate:.2%}/yr at every age and year"
        if self.mode == "1d":
            return f"1D age-varying scale from {self.source}"
        return (f"2D age x year scale from {self.source} "
                f"(years {self.years[0]}-{self.years[-1]}, held constant outside)")


def flat_scale(rate: float = 0.01) -> ImprovementScale:
    """Uniform improvement rate at every age, year, and sex."""
    return ImprovementScale(mode="flat", flat_rate=rate)


def _check_rates(values: np.ndarray, path: str) -> None:
    if np.isnan(values).any():
        raise ValueError(f"{path}: improvement table contains missing values.")
    if (np.abs(values) > _MAX_SANE_RATE).any():
        raise ValueError(
            f"{path}: found improvement rates with |rate| > {_MAX_SANE_RATE}. "
            "Rates must be decimals (0.01 = 1% per year), not percents."
        )


def _reindex_ages(series_by_age: pd.Series, max_age: int = 130) -> np.ndarray:
    """Return an array indexed 0..max_age, extending edge ages outward."""
    s = series_by_age.sort_index()
    full = s.reindex(range(0, max_age + 1))
    full = full.ffill().bfill()  # hold youngest/oldest rate beyond table edges
    return full.to_numpy(dtype=float)


def load_improvement_csv(path: str) -> ImprovementScale:
    """
    Load a 1D or 2D improvement scale from CSV, auto-detecting the schema
    from the header (any 4-digit column name => 2D wide format).
    """
    df = pd.read_csv(path)
    df.columns = [str(c).strip().lower() for c in df.columns]

    if "age" not in df.columns:
        raise ValueError(f"{path}: improvement CSV must have an 'age' column.")

    year_cols = [c for c in df.columns if _YEAR_RE.match(c)]

    if year_cols:
        return _load_2d(df, year_cols, path)
    return _load_1d(df, path)


def _load_1d(df: pd.DataFrame, path: str) -> ImprovementScale:
    if "improvement" in df.columns:
        male = female = df.set_index("age")["improvement"]
    elif "male_improvement" in df.columns and "female_improvement" in df.columns:
        idx = df.set_index("age")
        male, female = idx["male_improvement"], idx["female_improvement"]
    else:
        raise ValueError(
            f"{path}: 1D improvement CSV needs either an 'improvement' column "
            "or both 'male_improvement' and 'female_improvement'."
        )
    rates = {
        "male": _reindex_ages(male),
        "female": _reindex_ages(female),
    }
    for arr in rates.values():
        _check_rates(arr, path)
    return ImprovementScale(mode="1d", rates_1d=rates, source=path)


def _load_2d(df: pd.DataFrame, year_cols: list, path: str) -> ImprovementScale:
    years = sorted(int(c) for c in year_cols)
    if [str(y) for y in years] != [c for c in df.columns if _YEAR_RE.match(c)]:
        # Column order isn't required to be sorted; we sort internally.
        pass

    def frame_to_array(frame: pd.DataFrame) -> np.ndarray:
        by_age = frame.set_index("age")[[str(y) for y in years]].sort_index()
        by_age = by_age[~by_age.index.duplicated(keep="first")]
        full = by_age.reindex(range(0, 131)).ffill().bfill()
        arr = full.to_numpy(dtype=float)
        _check_rates(arr, path)
        return arr

    rates = {}
    if "sex" in df.columns:
        df["sex"] = df["sex"].astype(str).str.strip().str.lower()
        for sex in ("male", "female"):
            sub = df[df["sex"] == sex]
            if not sub.empty:
                rates[sex] = frame_to_array(sub)
    else:
        arr = frame_to_array(df)
        rates = {"male": arr, "female": arr}

    if not rates:
        raise ValueError(f"{path}: 2D CSV has a 'sex' column but no 'male'/'female' rows.")

    return ImprovementScale(mode="2d", rates_2d=rates, years=years, source=path)


def load_improvement_scale(path: Optional[str] = None,
                           flat_rate: float = 0.01) -> ImprovementScale:
    """CSV scale if a path is given, otherwise a flat scale at flat_rate."""
    if path is not None:
        return load_improvement_csv(path)
    return flat_scale(flat_rate)


def build_cohort_table(ssa_table: pd.DataFrame, starting_age: int,
                       scale: ImprovementScale, base_year: int,
                       max_age: int = 119) -> pd.DataFrame:
    """
    Bake generational mortality improvement into a copy of the life table for
    a cohort aged `starting_age` in `base_year`.

    The cohort reaches age x in calendar year base_year + (x - starting_age),
    so qx at age x gets (x - starting_age) years of compounded improvement.
    The terminal age keeps qx = 1.0.
    """
    table = ssa_table.copy()
    for sex in ("male", "female"):
        col = f"{sex}_qx"
        if col not in table.columns:
            continue
        for x in range(starting_age + 1, max_age + 1):
            factor = scale.cumulative_factor(x, base_year, x - starting_age, sex)
            table.loc[x, col] = min(1.0, float(table.loc[x, col]) * factor)
        table.loc[max_age, col] = 1.0
    return table


# ---------------------------------------------------------------------------
# Templates for users to fill out
# ---------------------------------------------------------------------------

def write_template_1d(path: str = "improvement_template_1d.csv",
                      rate: float = 0.01, max_age: int = 119) -> str:
    """Write a 1D template: one improvement rate per age, pre-filled."""
    df = pd.DataFrame({
        "age": range(0, max_age + 1),
        "male_improvement": rate,
        "female_improvement": rate,
    })
    df.to_csv(path, index=False)
    return path


def write_template_2d(path: str = "improvement_template_2d.csv",
                      rate: float = 0.01, max_age: int = 119,
                      first_year: int = 2023, n_years: int = 50) -> str:
    """
    Write a 2D template: rows = ages, columns = calendar years, pre-filled.
    Years after the last column hold that column's rate forever, so 50
    columns is plenty even for an 80-year projection.
    """
    data = {"age": list(range(0, max_age + 1))}
    for y in range(first_year, first_year + n_years):
        data[str(y)] = [rate] * (max_age + 1)
    pd.DataFrame(data).to_csv(path, index=False)
    return path
