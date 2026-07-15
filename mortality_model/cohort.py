"""
Core cohort death distribution model.

For a cohort at a given (age, sex):
  1. Compute d(t): unconditional probability of dying in year t (sums to 1.0)
  2. Apply COVID pullforward via f(t): remove fraction f(t) of d(t) deaths
  3. Normalize residual to get g(t): death distribution for COVID survivors
  4. Derive effective qx and mortality multiples vs the original SSA table
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from config import ModelConfig
from mortality_model.ssa_table import get_qx
from mortality_model.pullforward import build_pullforward_vector


def compute_death_distribution(
    starting_age: int,
    sex: str,
    ssa_table: pd.DataFrame,
    max_age: int = 119,
) -> np.ndarray:
    """
    Compute d(t): unconditional annual death probabilities for a cohort.

    d[i] = P(die in year i+1 | alive at starting_age), where i is 0-indexed.
    Sum of d equals 1.0 (everyone eventually dies).

    Parameters
    ----------
    starting_age : int
        Age of the cohort at the start of the projection (end of COVID).
    sex : str
        'male' or 'female'
    ssa_table : pd.DataFrame
    max_age : int

    Returns
    -------
    np.ndarray of shape (max_age - starting_age + 1,)
    """
    n_years = max_age - starting_age + 1
    d = np.zeros(n_years)

    survival = 1.0
    for i in range(n_years):
        age = starting_age + i
        qx = get_qx(ssa_table, age, sex)
        d[i] = survival * qx
        survival *= 1.0 - qx

    total = d.sum()
    if total > 0:
        d /= total  # normalize for numerical precision
    return d


def apply_pullforward(
    d: np.ndarray,
    starting_age: int,
    sex: str,
    cfg: ModelConfig,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Apply COVID pullforward to the death distribution.

    Returns
    -------
    d_remaining : np.ndarray
        Unconditional death probabilities after removing pulled-forward deaths.
        These are NOT normalized — their sum is the fraction of the cohort
        that survived COVID.
    f : np.ndarray
        Pullforward fractions for each year (same shape as d).
    """
    f = np.array(build_pullforward_vector(len(d), starting_age, cfg.pullforward))
    d_remaining = d * (1.0 - f)
    return d_remaining, f


def compute_survivor_distribution(
    d_remaining: np.ndarray,
) -> tuple[np.ndarray, float]:
    """
    Compute g(t): conditional death distribution for COVID survivors.

    Parameters
    ----------
    d_remaining : np.ndarray
        Unconditional remaining death probabilities (from apply_pullforward).

    Returns
    -------
    g : np.ndarray
        Conditional death distribution normalized so g.sum() = 1.0.
    survivor_fraction : float
        Fraction of the original cohort that survived COVID.
    """
    survivor_fraction = d_remaining.sum()
    if survivor_fraction <= 0:
        raise ValueError(
            "No survivors: the pullforward removes 100% of the cohort's "
            "remaining deaths (e.g. step shape with peak=1.0 and a horizon "
            "covering the whole remaining lifetime). Post-COVID LE is "
            "undefined for an empty cohort -- lower the peak fraction or "
            "shorten the horizon."
        )
    g = d_remaining / survivor_fraction
    return g, survivor_fraction


def compute_effective_qx(
    g: np.ndarray,
    starting_age: int,
    ssa_table: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute effective qx for COVID survivors at each future age.

    effective_qx(t) = g(t) / S_g(t-1)
    where S_g(t) = 1 - cumsum(g)[t] is the survivor's survival function.

    Returns
    -------
    effective_qx : np.ndarray
    ages : np.ndarray
        Corresponding ages (starting_age, starting_age+1, ...).
    """
    n = len(g)
    G = np.cumsum(g)  # CDF of g

    # Survival function for survivors: P(T > t) using 0-indexed t
    # S_g[0] = 1.0 (alive at t=0, i.e. before any projection year)
    # S_g[i] = 1 - G[i-1] for i >= 1  (survived through year i)
    S_g = np.ones(n)
    S_g[1:] = 1.0 - G[:-1]

    # Effective qx: P(die in year t | survived to year t-1) for survivors
    effective_qx = np.where(S_g > 0, g / S_g, np.nan)

    ages = np.arange(starting_age, starting_age + n)
    return effective_qx, ages


def compute_mortality_multiples(
    effective_qx: np.ndarray,
    starting_age: int,
    sex: str,
    ssa_table: pd.DataFrame,
) -> np.ndarray:
    """
    Compute mortality multiple vs SSA table: effective_qx / original_qx.

    Values < 1.0 mean survivors have lower mortality than the SSA table.

    Returns
    -------
    np.ndarray of same length as effective_qx.
    """
    n = len(effective_qx)
    multiples = np.full(n, np.nan)

    for i in range(n):
        age = starting_age + i
        orig_qx = get_qx(ssa_table, age, sex)
        eff = effective_qx[i]
        if orig_qx > 0 and not np.isnan(eff):
            multiples[i] = eff / orig_qx

    return multiples
