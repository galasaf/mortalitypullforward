"""
SSA Period Life Table loader.

Provides qx (annual probability of death) by single year of age (0–119) and sex.
Source approximation: SSA 2019 Period Life Table (pre-COVID baseline).

If a custom CSV is provided (columns: age, male_qx, female_qx), it is used instead.
"""
import math
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Anchor points from SSA 2019 Period Life Table
# (ages not listed are log-linearly interpolated)
# ---------------------------------------------------------------------------

_MALE_ANCHORS = {
    0: 0.00576,
    1: 0.000374,
    10: 0.000100,
    20: 0.001163,
    30: 0.001576,
    40: 0.003097,
    50: 0.007254,
    60: 0.015990,
    65: 0.020980,
    70: 0.031590,
    75: 0.048220,
    80: 0.074780,
    85: 0.117210,
    90: 0.181680,
    95: 0.269820,
    100: 0.377500,
    105: 0.500000,
    110: 0.650000,
    115: 0.850000,
    119: 1.000000,
}

_FEMALE_ANCHORS = {
    0: 0.00491,
    1: 0.000307,
    10: 0.000090,
    20: 0.000490,
    30: 0.000830,
    40: 0.001790,
    50: 0.004360,
    60: 0.009730,
    65: 0.013610,
    70: 0.021460,
    75: 0.034320,
    80: 0.057270,
    85: 0.097010,
    90: 0.159880,
    95: 0.246960,
    100: 0.360760,
    105: 0.480000,
    110: 0.620000,
    115: 0.800000,
    119: 1.000000,
}


def _interpolate_qx(anchors: dict, max_age: int = 119) -> np.ndarray:
    """
    Log-linearly interpolate between anchor points to produce qx for all ages 0–max_age.
    Uses linear interpolation in log space (i.e., geometric interpolation).
    """
    ages = sorted(anchors.keys())
    qx = np.zeros(max_age + 1)

    for i in range(len(ages) - 1):
        a0, a1 = ages[i], ages[i + 1]
        q0, q1 = anchors[a0], anchors[a1]

        log_q0 = math.log(q0)
        log_q1 = math.log(q1)

        for age in range(a0, a1 + 1):
            t = (age - a0) / (a1 - a0)
            qx[age] = math.exp(log_q0 + t * (log_q1 - log_q0))

    # Ensure the final age has qx = 1
    qx[max_age] = 1.0
    return qx


def _build_embedded_table(max_age: int = 119) -> pd.DataFrame:
    """Build the embedded SSA 2019 approximation as a DataFrame."""
    male_qx = _interpolate_qx(_MALE_ANCHORS, max_age)
    female_qx = _interpolate_qx(_FEMALE_ANCHORS, max_age)

    df = pd.DataFrame(
        {"male_qx": male_qx, "female_qx": female_qx},
        index=pd.RangeIndex(max_age + 1, name="age"),
    )
    return df


def load_ssa_table(path: Optional[str] = None, max_age: int = 119) -> pd.DataFrame:
    """
    Load SSA life table.

    Parameters
    ----------
    path : str or None
        Path to CSV with columns: age, male_qx, female_qx.
        If None, the built-in SSA 2019 approximation is used.
    max_age : int
        Maximum age to include. Defaults to 119.

    Returns
    -------
    pd.DataFrame
        Index = age (0 to max_age), columns = [male_qx, female_qx].
        qx values are clipped to [0, 1]; age max_age always has qx = 1.0.
    """
    if path is not None:
        df = pd.read_csv(path, index_col="age")
        df = df.reindex(range(max_age + 1))
        # Forward-fill any missing ages (shouldn't normally happen)
        df = df.ffill()
    else:
        df = _build_embedded_table(max_age)

    # Clip to valid probability range
    df["male_qx"] = df["male_qx"].clip(0.0, 1.0)
    df["female_qx"] = df["female_qx"].clip(0.0, 1.0)

    # Enforce terminal age
    df.loc[max_age, "male_qx"] = 1.0
    df.loc[max_age, "female_qx"] = 1.0

    return df


def get_qx(table: pd.DataFrame, age: int, sex: str) -> float:
    """Return qx for a given age and sex from the life table."""
    col = f"{sex}_qx"
    age = min(age, table.index.max())
    return float(table.loc[age, col])
