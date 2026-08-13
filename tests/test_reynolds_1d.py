"""1-D finite-width (short-bearing) Reynolds solver, validated against Ocvirk."""

import math

import pytest

from mochi.ocvirk_bearing import short_bearing_force, static_load_capacity_n
from mochi.reynolds_1d import solve_short_bearing_1d


@pytest.mark.parametrize("eps", [0.1, 0.3, 0.5, 0.7, 0.9])
def test_static_force_matches_ocvirk_closed_form(eps) -> None:
    """Pure rotation: the numeric F_e/F_t reproduce the Section 4.9 Ocvirk force."""

    omega = 95.0
    numeric = solve_short_bearing_1d(eps, omega)
    closed = short_bearing_force(eps, omega)
    scale = closed.magnitude_n
    assert abs(numeric.radial_force_n - closed.radial_force_n) / scale < 1.0e-3
    assert abs(numeric.tangential_force_n - closed.tangential_force_n) / scale < 1.0e-3


@pytest.mark.parametrize("eps", [0.2, 0.6, 0.85])
def test_magnitude_matches_load_capacity(eps) -> None:
    omega = 90.0
    numeric = solve_short_bearing_1d(eps, omega)
    assert numeric.magnitude_n == pytest.approx(static_load_capacity_n(eps, omega), rel=1.0e-3)


def test_radial_component_supports_the_load() -> None:
    """Under pure rotation F_e > 0 (restoring) and F_t < 0 (rotation sense)."""

    numeric = solve_short_bearing_1d(0.5, 100.0)
    assert numeric.radial_force_n > 0.0
    assert numeric.tangential_force_n < 0.0


def test_centred_journal_carries_no_load() -> None:
    numeric = solve_short_bearing_1d(0.0, 100.0)
    assert numeric.magnitude_n == pytest.approx(0.0, abs=1.0e-9)


def test_force_scales_linearly_with_speed() -> None:
    """Pure-rotation force is linear in the entrainment speed."""

    one = solve_short_bearing_1d(0.6, 50.0)
    two = solve_short_bearing_1d(0.6, 100.0)
    assert two.magnitude_n == pytest.approx(2.0 * one.magnitude_n, rel=1.0e-6)


def test_finer_grid_reduces_the_ocvirk_error() -> None:
    """Refining the circumferential grid moves the numeric force toward Ocvirk."""

    eps, omega = 0.7, 95.0
    closed = short_bearing_force(eps, omega).magnitude_n
    coarse = abs(solve_short_bearing_1d(eps, omega, n_circumferential=48).magnitude_n - closed)
    fine = abs(solve_short_bearing_1d(eps, omega, n_circumferential=720).magnitude_n - closed)
    assert fine < coarse


def test_squeeze_term_adds_radial_stiffness() -> None:
    """A closing film (eps_dot > 0) raises the radial force above pure rotation."""

    rotation = solve_short_bearing_1d(0.5, 100.0)
    squeezing = solve_short_bearing_1d(0.5, 100.0, eccentricity_rate_per_s=5.0)
    assert squeezing.radial_force_n > rotation.radial_force_n


def test_invalid_inputs_raise() -> None:
    with pytest.raises(ValueError):
        solve_short_bearing_1d(1.0, 100.0)  # eps must be < 1
    with pytest.raises(ValueError):
        solve_short_bearing_1d(-0.1, 100.0)
    with pytest.raises(ValueError):
        solve_short_bearing_1d(0.5, 100.0, viscosity_pa_s=0.0)
    with pytest.raises(ValueError):
        solve_short_bearing_1d(0.5, 100.0, n_axial=40)  # must be odd
    with pytest.raises(ValueError):
        solve_short_bearing_1d(0.5, 100.0, n_circumferential=4)


def test_matches_analytic_axial_parabola_peak() -> None:
    """Sanity: the closed form and numeric agree on the load angle (attitude)."""

    eps, omega = 0.5, 95.0
    numeric = solve_short_bearing_1d(eps, omega)
    closed = short_bearing_force(eps, omega)
    numeric_attitude = math.atan2(numeric.tangential_force_n, numeric.radial_force_n)
    assert numeric_attitude == pytest.approx(closed.attitude_angle_rad, abs=1.0e-3)
