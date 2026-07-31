"""Rotor lateral equation of motion and its time integration (D2 dynamics)."""

import math

import pytest
from scipy.integrate import solve_ivp

from mochi.chambers import build_cycle_trace
from mochi.journal_bearing import JOURNAL_CLEARANCE_M
from mochi.kinematics import RotaryGeometry
from mochi.ocvirk_bearing import (
    crank_pin_entrainment_speed_rad_s,
    equilibrium_eccentricity_ratio,
    short_bearing_force,
)
from mochi.rotor_dynamics import ROTOR_MASS_KG, _film_force_on_rotor, integrate_rotor_orbit
from mochi.true_gas_force import true_gas_load

TRACE_SAMPLES = 72
TRACE_STEP_M = 2.0e-4


@pytest.fixture(scope="module")
def trace():
    return build_cycle_trace(
        RotaryGeometry.default(), samples=TRACE_SAMPLES, grid_step_m=TRACE_STEP_M
    )


@pytest.fixture(scope="module")
def orbit(trace):
    return integrate_rotor_orbit(
        RotaryGeometry.default(),
        revolutions=3,
        samples=120,
        grid_samples=180,
        trace=trace,
        rtol=1.0e-6,
    )


@pytest.mark.parametrize("theta_deg", [150.0, 200.0, 226.0])
def test_static_film_force_cancels_the_gas_load(theta_deg, trace) -> None:
    """At the quasi-static equilibrium (e_dot = 0) the film reaction equals -F_gas."""

    geometry = RotaryGeometry.default()
    clearance = JOURNAL_CLEARANCE_M
    theta = math.radians(theta_deg)
    gas = true_gas_load(geometry, theta, trace=trace).rotor_force_n
    load = math.hypot(*gas)
    entrainment = crank_pin_entrainment_speed_rad_s(geometry, theta)
    eps = equilibrium_eccentricity_ratio(load, entrainment, clearance_m=clearance)
    psi = short_bearing_force(eps, entrainment, clearance_m=clearance).attitude_angle_rad
    direction = math.atan2(gas[1], gas[0]) - psi
    film = _film_force_on_rotor(
        eps * clearance * math.cos(direction),
        eps * clearance * math.sin(direction),
        0.0,
        0.0,
        entrainment,
        clearance,
    )
    residual = math.hypot(gas[0] + film[0], gas[1] + film[1])
    assert residual < 1.0e-6 * load


def test_frozen_load_settles_to_quasi_static(trace) -> None:
    """Reduction check: under a constant load the orbit relaxes to the quasi-static ε."""

    geometry = RotaryGeometry.default()
    clearance = JOURNAL_CLEARANCE_M
    mass = ROTOR_MASS_KG
    theta = math.radians(200.0)
    gas = true_gas_load(geometry, theta, trace=trace).rotor_force_n
    entrainment = crank_pin_entrainment_speed_rad_s(geometry, theta)
    eps_qs = equilibrium_eccentricity_ratio(math.hypot(*gas), entrainment, clearance_m=clearance)

    def rhs(_t, state):
        e_x, e_y, v_x, v_y = state
        film = _film_force_on_rotor(e_x, e_y, v_x, v_y, entrainment, clearance)
        return [v_x, v_y, (gas[0] + film[0]) / mass, (gas[1] + film[1]) / mass]

    solution = solve_ivp(
        rhs, (0.0, 0.5), [clearance * 0.2, 0.0, 0.0, 0.0], method="BDF", rtol=1e-8, atol=1e-13
    )
    eps_dyn = math.hypot(solution.y[0, -1], solution.y[1, -1]) / clearance
    assert eps_dyn == pytest.approx(eps_qs, abs=1.0e-3)


def test_orbit_is_bounded_and_squeeze_attenuates_the_peak(orbit) -> None:
    """The dynamic peak stays below the quasi-static estimate (squeeze-film lag)."""

    assert 0.0 < orbit.peak_eccentricity_ratio < 1.0
    assert orbit.peak_eccentricity_ratio < orbit.quasi_static_peak_eccentricity_ratio
    # a larger minimum film than the quasi-static (over-)estimate of ~4.4 um
    assert orbit.minimum_film_thickness_m > 5.0e-6


def test_orbit_is_periodic(orbit) -> None:
    """The steady revolution closes on itself (same crank angle, same eccentricity)."""

    tolerance = 0.02 * orbit.clearance_m
    assert orbit.eccentricity_x_m[0] == pytest.approx(orbit.eccentricity_x_m[-1], abs=tolerance)
    assert orbit.eccentricity_y_m[0] == pytest.approx(orbit.eccentricity_y_m[-1], abs=tolerance)


def test_dynamic_journal_friction_vs_quasi_static(orbit) -> None:
    """The dynamic journal friction exceeds concentric Petroff and tracks quasi-static."""

    assert orbit.petroff_journal_friction_w > 0.0
    # eccentric (dynamic) friction is above the concentric Petroff estimate
    assert orbit.dynamic_journal_friction_w > orbit.petroff_journal_friction_w
    # the attenuated dynamic eccentricity keeps it close to (slightly under) quasi-static
    assert orbit.dynamic_journal_friction_w == pytest.approx(
        orbit.quasi_static_journal_friction_w, rel=0.1
    )


def test_centrifugal_load_is_m_omega2_e(orbit) -> None:
    geometry = RotaryGeometry.default()
    expected = ROTOR_MASS_KG * geometry.angular_speed_rad_s**2 * geometry.eccentricity_m
    assert orbit.centrifugal_load_n == pytest.approx(expected, rel=1.0e-9)


def test_zero_eccentricity_film_force_is_zero() -> None:
    assert _film_force_on_rotor(0.0, 0.0, 0.0, 0.0, 95.0, JOURNAL_CLEARANCE_M) == (0.0, 0.0)


def test_invalid_inputs_raise() -> None:
    geometry = RotaryGeometry.default()
    with pytest.raises(ValueError):
        integrate_rotor_orbit(geometry, mass_kg=0.0)
    with pytest.raises(ValueError):
        integrate_rotor_orbit(geometry, revolutions=0)
