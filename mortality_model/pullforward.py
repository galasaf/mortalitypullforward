"""
COVID pullforward distribution.

f(t) = fraction of people who *would* have died in projection year t
       that instead died during COVID.

Convention: t is 1-indexed (t=1 = first year of projection from end of COVID).
"""
import math
from config import PullforwardConfig


def get_grade_out_years(age: int, cfg: PullforwardConfig) -> int:
    """
    Return the linear grade-out period (years) for a given starting age.
    Uses age-band overrides if cfg.age_varying is True; otherwise returns default.
    """
    if not cfg.age_varying:
        return cfg.default_grade_out_years

    for (min_age, max_age), years in cfg.age_grade_out.items():
        if min_age <= age < max_age:
            return years

    return cfg.default_grade_out_years


def compute_pullforward_fraction(t: int, age: int, cfg: PullforwardConfig) -> float:
    """
    Return f(t): fraction pulled forward for projection year t.

    Parameters
    ----------
    t : int
        Projection year (1-indexed). t=1 = highest severity.
    age : int
        Starting age of the cohort (used for age-varying grade-out).
    cfg : PullforwardConfig

    Returns
    -------
    float in [0, 1]
    """
    if t <= 0:
        return 0.0

    if cfg.shape == "linear":
        grade_out = get_grade_out_years(age, cfg)
        if grade_out <= 0:
            return 0.0  # grade-out 0 = no pullforward at all
        shape_val = max(0.0, 1.0 - (t - 1) / grade_out)
    elif cfg.shape == "step":
        # Constant box: full peak effect for every year inside the horizon.
        grade_out = get_grade_out_years(age, cfg)
        shape_val = 1.0 if t <= grade_out else 0.0
    elif cfg.shape == "exponential":
        shape_val = math.exp(-cfg.exponential_decay_rate * (t - 1))
    else:
        raise ValueError(f"Unknown pullforward shape: {cfg.shape!r}")

    return cfg.peak_fraction * shape_val


def build_pullforward_vector(n_years: int, age: int, cfg: PullforwardConfig) -> list:
    """
    Return a list f[0..n_years-1] where f[i] = pullforward fraction for year (i+1).
    """
    return [compute_pullforward_fraction(t + 1, age, cfg) for t in range(n_years)]
