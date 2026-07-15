"""
All configurable parameters for the mortality pullforward model.
Edit this file to run different scenarios, or use the named presets at the bottom.
"""
from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple


@dataclass
class PullforwardConfig:
    """
    Controls the COVID pullforward distribution f(t):
    the fraction of 'would-be year-t deaths' that instead died during COVID.

    f(t) = peak_fraction * shape_function(t)

    peak_fraction scales the entire curve. At t=1, shape_function=1.0, so
    peak_fraction is the fraction of 'year-1 deaths' that were pulled into COVID.
    Realistic range: 0.3–0.8. Higher = stronger harvesting effect.
    """

    # ---- Peak severity ----
    # Fraction of year-1 deaths pulled into COVID. Scales the entire f(t) curve.
    # 1.0 = 100% of people who would have died in year 1 died during COVID (aggressive)
    # 0.5 = 50% of year-1 deaths pulled forward (moderate)
    # 0.3 = 30% (conservative)
    peak_fraction: float = 1.0

    # ---- Shape ----
    # 'linear':      f(t) = peak_fraction * max(0, 1 - (t-1) / grade_out_years)
    # 'step':        f(t) = peak_fraction for t <= grade_out_years, 0 after
    #                [constant box: full effect for every year inside the horizon]
    # 'exponential': f(t) = peak_fraction * exp(-decay_rate * (t-1))  [fat tail, no hard cutoff]
    shape: str = "linear"

    # Linear grade-out: years over which pullforward grades from peak → 0
    # e.g. grade_out=8, peak=0.6: t=1→60%, t=5→22.5%, t=8→7.5%, t=9→0%
    # grade_out=0 disables the pullforward entirely (f(t)=0 for all t).
    default_grade_out_years: int = 10

    # Exponential decay rate (only used when shape='exponential')
    # Higher value = steeper / faster decay toward zero
    exponential_decay_rate: float = 0.3

    # ---- Age-varying grade-out ----
    # If True, older cohorts use a shorter grade-out. Rationale: older people have
    # shorter remaining LE, so COVID could only plausibly pull forward nearer-term deaths.
    age_varying: bool = True

    # Map: (min_age_inclusive, max_age_exclusive) → grade_out_years override.
    # Applied only when age_varying=True and shape is 'linear' or 'step'.
    age_grade_out: Dict[Tuple[int, int], int] = field(
        default_factory=lambda: {
            (0, 65):   10,
            (65, 75):   8,
            (75, 85):   6,
            (85, 200):  5,
        }
    )


@dataclass
class ModelConfig:
    """Top-level model configuration."""

    # ---- COVID timing ----
    covid_end_year: int = 2022  # Reference point: end of COVID / start of post-COVID LE

    # ---- SSA life table ----
    # Path to CSV with columns: age, male_qx, female_qx (ages 0–119).
    # Set to None to use the built-in embedded approximation of SSA 2019.
    ssa_table_path: Optional[str] = None

    # ---- Projection ----
    max_age: int = 119

    # ---- Analysis scope ----
    analysis_ages: list = field(
        default_factory=lambda: [40, 50, 55, 60, 65, 70, 75, 80]
    )
    analysis_sexes: list = field(default_factory=lambda: ["male", "female"])

    # ---- Mortality improvement ----
    # Future mortality improvement compounds from covid_end_year forward:
    # projection year t uses qx * (1 - rate)^(t-1) for flat scales, or the
    # equivalent cumulative product for age/year-varying CSV scales.
    improvement_enabled: bool = True

    # Flat annual improvement rate used when no CSV table is supplied.
    # 0.01 = 1% lower mortality per calendar year at every age (the default).
    improvement_rate: float = 0.01

    # Path to a 1D or 2D improvement CSV (schema auto-detected; see
    # mortality_model/improvement.py docstring, or generate a template with
    # `python sensitivity.py --improvement-template 1d|2d`). None = flat rate.
    improvement_table_path: Optional[str] = None

    # ---- Pullforward ----
    pullforward: PullforwardConfig = field(default_factory=PullforwardConfig)


# ---------------------------------------------------------------------------
# Named scenario presets — used by scenarios.py
#
# Dimensions being varied:
#   peak_fraction  : how severe the pullforward is at its peak
#   grade_out      : how many years out the effect reaches
#   age_varying    : whether older cohorts have steeper distributions
#   shape          : linear (hard cutoff) vs exponential (fat tail)
# ---------------------------------------------------------------------------

def default_config() -> ModelConfig:
    """
    Baseline: linear pullforward, age-varying grade-out, peak=1.0.
    Use as a reference point; peak=1.0 is aggressive.
    """
    return ModelConfig()


def short_harvest_config() -> ModelConfig:
    """
    Conservative: flu-harvesting literature baseline.
    COVID only pulled forward deaths within ~3 years; moderate peak.
    Interpretation: most COVID victims would have died within 3 years anyway,
    but not all of them.
    """
    pf = PullforwardConfig(
        shape="linear",
        peak_fraction=0.60,
        default_grade_out_years=3,
        age_varying=False,
    )
    return ModelConfig(pullforward=pf)


def moderate_base_config() -> ModelConfig:
    """
    Moderate: most defensible for COVID given its severity vs flu.
    7-year grade-out with age-varying bands and a realistic peak.
    """
    pf = PullforwardConfig(
        shape="linear",
        peak_fraction=0.65,
        default_grade_out_years=7,
        age_varying=True,
        age_grade_out={
            (0, 65):   10,
            (65, 75):   7,
            (75, 85):   5,
            (85, 200):  3,
        },
    )
    return ModelConfig(pullforward=pf)


def elderly_concentrated_config() -> ModelConfig:
    """
    Age-stratified: COVID's impact was concentrated on the frailest elderly.
    Very steep grade-out for 85+, mild for younger cohorts.
    Interpretation: COVID was essentially a near-term mortality accelerator
    only for those already very close to dying.
    """
    pf = PullforwardConfig(
        shape="linear",
        peak_fraction=0.70,
        age_varying=True,
        age_grade_out={
            (0, 65):   15,   # minimal pullforward for young (<65)
            (65, 75):   8,
            (75, 85):   5,
            (85, 200):  3,   # steep for 85+: nearly all near-term deaths harvested
        },
    )
    return ModelConfig(pullforward=pf)


def long_harvest_config() -> ModelConfig:
    """
    Aggressive: COVID pulled forward deaths that were a decade+ out.
    Represents the upper bound of plausible pullforward.
    """
    pf = PullforwardConfig(
        shape="linear",
        peak_fraction=0.50,
        default_grade_out_years=15,
        age_varying=True,
        age_grade_out={
            (0, 65):   15,
            (65, 75):  12,
            (75, 85):  10,
            (85, 200):  7,
        },
    )
    return ModelConfig(pullforward=pf)


def exponential_moderate_config() -> ModelConfig:
    """
    Exponential shape, no hard cutoff, moderate peak and decay.
    Fat tail: some pullforward extends indefinitely but diminishes quickly.
    Decay rate 0.4 means ~50% reduction every ~1.7 years.
    """
    pf = PullforwardConfig(
        shape="exponential",
        peak_fraction=0.70,
        exponential_decay_rate=0.4,
        age_varying=False,
    )
    return ModelConfig(pullforward=pf)


def exponential_long_tail_config() -> ModelConfig:
    """
    Exponential shape with a slow decay (fat tail).
    High peak, slow fade — represents a world where COVID's selection
    effect extended far into the future with no hard cutoff.
    Decay rate 0.2 means ~50% reduction every ~3.5 years.
    """
    pf = PullforwardConfig(
        shape="exponential",
        peak_fraction=0.80,
        exponential_decay_rate=0.2,
        age_varying=False,
    )
    return ModelConfig(pullforward=pf)


def sickest_only_config() -> ModelConfig:
    """
    Very short grade-out, high peak.
    Interpretation: COVID almost exclusively killed people who were within
    2 years of dying anyway — but captured nearly all of them.
    Similar to a sharp actuarial 'deaths pulled forward by <2 years' view.
    """
    pf = PullforwardConfig(
        shape="linear",
        peak_fraction=0.90,
        default_grade_out_years=2,
        age_varying=False,
    )
    return ModelConfig(pullforward=pf)


# ---------------------------------------------------------------------------
# Ad-hoc config builder — used by sensitivity.py to override parameters from
# the command line without editing this file.
# ---------------------------------------------------------------------------

def apply_overrides(
    cfg: "ModelConfig",
    peak: Optional[float] = None,
    grade_out: Optional[int] = None,
    shape: Optional[str] = None,
    decay: Optional[float] = None,
    age_varying: Optional[bool] = None,
    ages: Optional[list] = None,
    sexes: Optional[list] = None,
    improvement_rate: Optional[float] = None,
    improvement_table: Optional[str] = None,
    no_improvement: Optional[bool] = None,
) -> "ModelConfig":
    """
    Mutate `cfg` in place with any provided overrides and return it.

    Only non-None arguments take effect, so callers can override one knob at a
    time. Note on grade_out: a single grade-out value is flat across ages, so
    setting it disables age-varying bands unless `age_varying=True` is also
    passed explicitly.
    """
    pf = cfg.pullforward

    if peak is not None:
        pf.peak_fraction = peak
    if shape is not None:
        if shape not in ("linear", "step", "exponential"):
            raise ValueError(f"shape must be 'linear', 'step', or 'exponential', got {shape!r}")
        pf.shape = shape
    if decay is not None:
        pf.exponential_decay_rate = decay
    if grade_out is not None:
        pf.default_grade_out_years = grade_out
        # A flat grade-out only takes effect when age-varying bands are off.
        if age_varying is None:
            pf.age_varying = False
    if age_varying is not None:
        pf.age_varying = age_varying

    if ages is not None:
        cfg.analysis_ages = ages
    if sexes is not None:
        cfg.analysis_sexes = sexes

    if improvement_rate is not None:
        cfg.improvement_rate = improvement_rate
        cfg.improvement_enabled = True
    if improvement_table is not None:
        cfg.improvement_table_path = improvement_table
        cfg.improvement_enabled = True
    if no_improvement:
        cfg.improvement_enabled = False

    return cfg


def describe_config(cfg: "ModelConfig") -> str:
    """Return a one-block human-readable summary of the active assumptions."""
    pf = cfg.pullforward
    lines = [
        f"  Shape          : {pf.shape}",
        f"  Peak fraction  : {pf.peak_fraction:.0%}",
        f"  Age-varying    : {pf.age_varying}",
    ]
    if pf.shape in ("linear", "step"):
        if pf.age_varying:
            lines.append(f"  Grade-out bands: {pf.age_grade_out}")
        else:
            lines.append(f"  Grade-out years: {pf.default_grade_out_years}")
    else:
        lines.append(f"  Decay rate     : {pf.exponential_decay_rate}")
    lines.append(f"  Ages           : {cfg.analysis_ages}")
    lines.append(f"  Sexes          : {cfg.analysis_sexes}")
    lines.append(f"  Improvement    : {describe_improvement(cfg)}")
    return "\n".join(lines)


def describe_improvement(cfg: "ModelConfig") -> str:
    """One-line summary of the mortality improvement assumption."""
    if not cfg.improvement_enabled:
        return "OFF (static life table)"
    if cfg.improvement_table_path:
        return f"from CSV: {cfg.improvement_table_path}"
    return f"flat {cfg.improvement_rate:.2%}/yr at every age and year"


def uniform_benchmark_config() -> ModelConfig:
    """
    Clean theoretical benchmark: no age-variation, moderate assumptions.
    Useful as a baseline before adding age-stratification.
    """
    pf = PullforwardConfig(
        shape="linear",
        peak_fraction=0.50,
        default_grade_out_years=10,
        age_varying=False,
    )
    return ModelConfig(pullforward=pf)
