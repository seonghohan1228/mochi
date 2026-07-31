"""Short-bearing (Ocvirk) crank-pin journal force -- eccentricity-coupled film."""

import math

import pytest

from mochi.chambers import build_cycle_trace
from mochi.journal_bearing import (
    JOURNAL_CLEARANCE_M,
    JOURNAL_LENGTH_M,
    JOURNAL_RADIUS_M,
    LUBRICANT_VISCOSITY_PA_S,
    journal_relative_speed_rad_s,
    petroff_friction,
)
from mochi.kinematics import RotaryGeometry
from mochi.ocvirk_bearing import (
    crank_pin_entrainment_speed_rad_s,
    eccentric_friction_torque_nm,
    eccentricity_cycle,
    equilibrium_eccentricity_ratio,
    short_bearing_force,
    static_load_capacity_n,
)
from mochi.true_gas_force import peak_rotor_force_n

_PREFACTOR = (
    LUBRICANT_VISCOSITY_PA_S * JOURNAL_RADIUS_M * JOURNAL_LENGTH_M**3 / JOURNAL_CLEARANCE_M**2
)

# A coarse trace keeps the cycle sweeps cheap (the full 121-sample raster is slow);
# 2e-4 m grid and 72 samples match the other bearing/gas-force test modules.
TRACE_SAMPLES = 72
TRACE_STEP_M = 2.0e-4
SWEEP_SAMPLES = 90


@pytest.fixture(scope="module")
def trace():
    return build_cycle_trace(
        RotaryGeometry.default(), samples=TRACE_SAMPLES, grid_step_m=TRACE_STEP_M
    )


def _textbook_load(eps: float, omega: float) -> float:
    """Ocvirk pure-rotation load capacity |W| (independent closed form)."""

    return (
        _PREFACTOR
        * omega
        * eps
        * math.sqrt(16.0 * eps**2 + math.pi**2 * (1.0 - eps**2))
        / (2.0 * (1.0 - eps**2) ** 2)
    )


@pytest.mark.parametrize("eps", [0.1, 0.3, 0.5, 0.7, 0.9])
def test_static_load_matches_textbook_ocvirk_capacity(eps) -> None:
    """|W|(eps) equals the standard short-bearing load-capacity closed form."""

    omega = 95.0
    assert static_load_capacity_n(eps, omega) == pytest.approx(
        _textbook_load(eps, omega), rel=1e-12
    )


@pytest.mark.parametrize("eps", [0.2, 0.5, 0.8])
def test_attitude_angle_matches_short_bearing_relation(eps) -> None:
    """|tan psi| = pi sqrt(1-eps^2) / (4 eps) for pure rotation."""

    force = short_bearing_force(eps, 100.0)
    expected = math.pi * math.sqrt(1.0 - eps**2) / (4.0 * eps)
    assert abs(math.tan(force.attitude_angle_rad)) == pytest.approx(expected, rel=1e-9)


def test_force_components_match_closed_form() -> None:
    """Radial/tangential components equal the P, Q, S closed form (eps_dot = 0)."""

    eps, omega = 0.6, 80.0
    one_minus = 1.0 - eps**2
    q_int = -2.0 * eps / one_minus**2
    s_int = math.pi / (2.0 * one_minus**1.5)
    force = short_bearing_force(eps, omega)
    assert force.radial_force_n == pytest.approx(_PREFACTOR * (-eps * omega * q_int), rel=1e-12)
    assert force.tangential_force_n == pytest.approx(_PREFACTOR * (-eps * omega * s_int), rel=1e-12)
    assert force.magnitude_n == pytest.approx(
        math.hypot(force.radial_force_n, force.tangential_force_n)
    )


def test_radial_force_supports_the_load() -> None:
    """The line-of-centres component is restoring (positive) under rotation."""

    force = short_bearing_force(0.5, 100.0)
    assert force.radial_force_n > 0.0
    assert force.min_film_thickness_m == pytest.approx(JOURNAL_CLEARANCE_M * 0.5)


def test_squeeze_term_adds_radial_stiffness() -> None:
    """A closing film (eps_dot > 0) raises the radial force above pure rotation."""

    rotation_only = short_bearing_force(0.5, 100.0)
    squeezing = short_bearing_force(0.5, 100.0, eccentricity_rate_per_s=5.0)
    assert squeezing.radial_force_n > rotation_only.radial_force_n


def test_friction_reduces_to_petroff_at_zero_eccentricity() -> None:
    """As eps -> 0 the eccentric Couette torque recovers Petroff exactly (Section 4.7)."""

    geometry = RotaryGeometry.default()
    theta = math.radians(200.0)
    speed = journal_relative_speed_rad_s(geometry, theta)
    petroff = petroff_friction(geometry, theta).friction_torque_nm
    assert eccentric_friction_torque_nm(0.0, speed) == pytest.approx(petroff, rel=1e-12)


@pytest.mark.parametrize("eps", [0.3, 0.6, 0.85])
def test_eccentric_friction_is_petroff_over_sqrt(eps) -> None:
    """T_f(eps) = T_petroff / sqrt(1 - eps^2), so it grows with eccentricity."""

    speed = 180.0
    base = eccentric_friction_torque_nm(0.0, speed)
    assert eccentric_friction_torque_nm(eps, speed) == pytest.approx(
        base / math.sqrt(1.0 - eps**2), rel=1e-12
    )


def test_equilibrium_inverts_the_load_capacity() -> None:
    """equilibrium_eccentricity_ratio round-trips against static_load_capacity_n."""

    omega = 90.0
    eps = 0.55
    load = static_load_capacity_n(eps, omega)
    assert equilibrium_eccentricity_ratio(load, omega) == pytest.approx(eps, abs=1e-6)


def test_equilibrium_grows_with_load() -> None:
    omega = 90.0
    light = equilibrium_eccentricity_ratio(200.0, omega)
    heavy = equilibrium_eccentricity_ratio(2000.0, omega)
    assert 0.0 < light < heavy < 1.0


def test_zero_load_gives_centred_journal() -> None:
    assert equilibrium_eccentricity_ratio(0.0, 90.0) == 0.0


def test_unsupportable_load_clamps_to_contact() -> None:
    """A load beyond the film capacity at any eps < 1 clamps just below contact."""

    eps = equilibrium_eccentricity_ratio(1.0e9, 90.0)
    assert 0.999 < eps < 1.0


def test_entrainment_is_shaft_minus_half_relative() -> None:
    """Omega = omega - |omega_rel|/2 (mean surface speed)."""

    geometry = RotaryGeometry.default()
    theta = math.radians(140.0)
    omega = geometry.angular_speed_rad_s
    relative = journal_relative_speed_rad_s(geometry, theta)
    assert crank_pin_entrainment_speed_rad_s(geometry, theta) == pytest.approx(
        omega - 0.5 * abs(relative), rel=1e-12
    )


def test_cycle_peak_eccentricity_tracks_the_peak_reaction(trace) -> None:
    """The peak running eccentricity coincides with the true peak gas reaction (~2.5 kN)."""

    geometry = RotaryGeometry.default()
    cycle = eccentricity_cycle(geometry, samples=SWEEP_SAMPLES, trace=trace)
    assert 0.0 < cycle.peak_eccentricity_ratio < 1.0
    assert max(cycle.load_n) == pytest.approx(
        peak_rotor_force_n(geometry, samples=SWEEP_SAMPLES, trace=trace), rel=0.02
    )
    # The minimum film thickness is the clearance scaled by the peak eccentricity.
    assert cycle.minimum_film_thickness_m == pytest.approx(
        JOURNAL_CLEARANCE_M * (1.0 - cycle.peak_eccentricity_ratio), rel=1e-9
    )


def test_cycle_friction_exceeds_petroff(trace) -> None:
    """Running off-centre dissipates more than the concentric Petroff estimate."""

    cycle = eccentricity_cycle(RotaryGeometry.default(), samples=SWEEP_SAMPLES, trace=trace)
    assert cycle.mean_friction_power_w > cycle.petroff_mean_friction_power_w > 0.0


def test_invalid_inputs_raise() -> None:
    with pytest.raises(ValueError):
        short_bearing_force(1.0, 100.0)  # eps must be < 1
    with pytest.raises(ValueError):
        short_bearing_force(-0.1, 100.0)
    with pytest.raises(ValueError):
        short_bearing_force(0.5, 100.0, viscosity_pa_s=0.0)
    with pytest.raises(ValueError):
        eccentric_friction_torque_nm(1.5, 100.0)
    with pytest.raises(ValueError):
        equilibrium_eccentricity_ratio(-5.0, 100.0)
    with pytest.raises(ValueError):
        # The samples guard fires before any trace is built, so this stays cheap.
        eccentricity_cycle(RotaryGeometry.default(), samples=4)
