"""Rotor-cylinder sealing-line sliding kinematics and boundary friction."""

import math

import numpy as np
import pytest

from mochi.chambers import build_cycle_trace
from mochi.kinematics import RotaryGeometry
from mochi.rotor_cylinder import (
    ROTOR_CYLINDER_FRICTION_COEFF,
    contact_normal_force_n,
    contact_sliding_speed_m_s,
    cycle_mean_abs_sliding_speed_m_s,
    ehl_film_thickness_m,
    hertz_line_contact_force_n,
    integrate_sealing_contact,
    mixed_ehl_friction_power_w,
    rotor_cylinder_friction_power_w,
)


@pytest.fixture(scope="module")
def geometry():
    return RotaryGeometry.default()


@pytest.fixture(scope="module")
def trace(geometry):
    return build_cycle_trace(geometry, samples=48, grid_step_m=3.0e-4)


@pytest.fixture(scope="module")
def seal_orbit(geometry, trace):
    return integrate_sealing_contact(
        geometry, revolutions=2, samples=120, grid_samples=90, trace=trace
    )


def test_swing_rotor_slides_not_rolls(geometry):
    """Net-zero rotation over a cycle => the sealing line rubs (mean |v| ~ orbit speed)."""

    mean_abs = cycle_mean_abs_sliding_speed_m_s(geometry)
    orbit = geometry.eccentricity_m * geometry.angular_speed_rad_s  # e*omega ~ 0.85 m/s
    # A true rolling piston would give mean|v| ~ 0; the swing rotor is O(e*omega).
    assert mean_abs > 0.5 * orbit
    assert mean_abs == pytest.approx(0.95, abs=0.1)


def test_sliding_speed_is_orbit_plus_swing(geometry):
    """v_slide = -e*omega + omega*(dphi/dtheta)*R_r.

    The orientation swing peaks at theta = 90 deg / 270 deg, so d(phi)/d(theta) = 0
    there and the sliding is pure orbital, v = -e*omega. (At TDC/BDC the swing rate is
    at its maximum instead, so the seal slides fastest there, up to ~2 m/s.)
    """

    omega = geometry.angular_speed_rad_s
    orbit = -geometry.eccentricity_m * omega
    assert contact_sliding_speed_m_s(geometry, math.pi / 2) == pytest.approx(orbit, rel=1e-3)
    assert contact_sliding_speed_m_s(geometry, 3 * math.pi / 2) == pytest.approx(orbit, rel=1e-3)


def test_friction_scales_with_load_and_coefficient(geometry):
    """Boundary friction is linear in N_c and in mu."""

    base = rotor_cylinder_friction_power_w(geometry, 100.0, friction_coefficient=0.1)
    assert rotor_cylinder_friction_power_w(
        geometry, 200.0, friction_coefficient=0.1
    ) == pytest.approx(2.0 * base, rel=1e-9)
    assert rotor_cylinder_friction_power_w(
        geometry, 100.0, friction_coefficient=0.2
    ) == pytest.approx(2.0 * base, rel=1e-9)


def test_friction_equals_mu_load_meanspeed(geometry):
    """Constant-load friction equals mu * N_c * <|v_slide|>."""

    mu, load = 0.1, 100.0
    expected = mu * load * cycle_mean_abs_sliding_speed_m_s(geometry)
    assert rotor_cylinder_friction_power_w(
        geometry, load, friction_coefficient=mu
    ) == pytest.approx(expected, rel=1e-9)


def test_callable_load_profile(geometry):
    """A theta-dependent N_c(theta) is accepted (for a future N_c estimate)."""

    power = rotor_cylinder_friction_power_w(geometry, lambda th: 100.0 * (1.0 + 0.5 * math.cos(th)))
    assert power > 0.0


def test_default_coefficient_is_boundary_value():
    assert 0.04 <= ROTOR_CYLINDER_FRICTION_COEFF <= 0.15


def test_invalid_inputs_raise(geometry):
    with pytest.raises(ValueError):
        rotor_cylinder_friction_power_w(geometry, -10.0)
    with pytest.raises(ValueError):
        rotor_cylinder_friction_power_w(geometry, 100.0, friction_coefficient=-0.1)
    with pytest.raises(ValueError):
        contact_sliding_speed_m_s(geometry, float("nan"))


# ---------------------------------------------------------------------------
# Contact normal force N_c and the boundary-vs-mixed-EHL friction comparison
# ---------------------------------------------------------------------------


def test_contact_normal_force_has_centrifugal_floor(geometry, trace):
    """N_c is non-negative, O(tens of N) (journal takes the gas load), zero in seal-over."""

    m_omega2_e = 0.275 * geometry.angular_speed_rad_s**2 * geometry.eccentricity_m  # ~44 N
    forces = [
        contact_normal_force_n(geometry, 2.0 * math.pi * i / 120, trace=trace) for i in range(120)
    ]
    assert all(f >= 0.0 for f in forces)
    assert 20.0 < np.mean(forces) < 150.0  # residual after the journal, order tens of N
    assert max(forces) > m_omega2_e  # gas radial adds to the centrifugal floor
    # theta = 0 lies inside the seal-over window -> no defined contact.
    assert contact_normal_force_n(geometry, 0.0, trace=trace) == 0.0


def test_mixed_ehl_is_below_pure_boundary(geometry, trace):
    """A partial EHL film (conformal contact) carries load -> less friction than boundary."""

    load = lambda th: contact_normal_force_n(geometry, th, trace=trace)  # noqa: E731
    boundary = rotor_cylinder_friction_power_w(
        geometry, load, friction_coefficient=0.10, samples=48
    )
    mixed = mixed_ehl_friction_power_w(
        geometry, load, roughness_m=0.2e-6, boundary_coefficient=0.10, samples=48
    )
    assert 0.0 < mixed < boundary


def test_mixed_ehl_rougher_surface_approaches_boundary(geometry, trace):
    """Rougher surfaces (Lambda -> 0) push the mixed loss toward the pure-boundary bound."""

    load = lambda th: contact_normal_force_n(geometry, th, trace=trace)  # noqa: E731
    smooth = mixed_ehl_friction_power_w(geometry, load, roughness_m=0.05e-6, samples=48)
    rough = mixed_ehl_friction_power_w(geometry, load, roughness_m=1.0e-6, samples=48)
    assert smooth < rough


def test_mixed_ehl_invalid_inputs_raise(geometry):
    with pytest.raises(ValueError):
        mixed_ehl_friction_power_w(geometry, 100.0, roughness_m=0.0)
    with pytest.raises(ValueError):
        mixed_ehl_friction_power_w(geometry, 100.0, viscosity_pa_s=-1.0)


# ---------------------------------------------------------------------------
# Self-consistent sealing-contact rung (compliant Hertz contact in the EOM)
# ---------------------------------------------------------------------------


def test_hertz_line_contact_load_deflection(geometry):
    """Palmgren steel line contact: zero when apart, monotone, ~100 N at 0.22 um."""

    assert hertz_line_contact_force_n(-1.0e-6, 0.021) == 0.0
    steep = hertz_line_contact_force_n(0.5e-6, 0.021)
    shallow = hertz_line_contact_force_n(0.2e-6, 0.021)
    assert steep > shallow > 0.0
    assert hertz_line_contact_force_n(0.217e-6, 0.021) == pytest.approx(100.0, rel=0.1)


def test_self_consistent_collapses_penetration(seal_orbit):
    """The ~6 um free-orbit penetration collapses to the physical Hertz deflection."""

    assert 0.0 < seal_orbit.max_penetration_m < 3.0e-6  # was ~6 um without the reaction
    assert seal_orbit.mean_normal_force_n > 0.0
    assert seal_orbit.peak_normal_force_n > seal_orbit.mean_normal_force_n


def test_self_consistent_friction_boundary_vs_mixed(seal_orbit):
    """Both closures are positive and comparable.

    In the mixed / full-film regime the hydrodynamic film unloads the asperities, so
    mixed-EHL sits *below* boundary; in the boundary regime (Lambda < 1, the Ra 0.3 um
    design finish) the asperities carry the whole load and mixed-EHL converges to
    boundary plus a small viscous-shear term. So they agree to within ~15 %.
    """

    assert seal_orbit.boundary_friction_power_w > 0.0
    assert seal_orbit.mixed_ehl_friction_power_w > 0.0
    assert seal_orbit.mixed_ehl_friction_power_w == pytest.approx(
        seal_orbit.boundary_friction_power_w, rel=0.15
    )


def test_ehl_film_puts_the_seal_in_the_mixed_regime(seal_orbit):
    """The hard-EHL film is thin (Lambda ~ 1) -- the seal does NOT float, it runs mixed."""

    assert 0.0 < seal_orbit.mean_ehl_film_thickness_m < 1.0e-6  # sub-micron, not a thick film
    assert seal_orbit.mean_film_parameter < 3.0  # mixed / boundary, not full-film


def test_ehl_film_thickness_thickens_with_speed_thins_with_load(geometry):
    """Dowson-Higginson: h_min rises with entrainment, falls with load; inf when unloaded."""

    fast = ehl_film_thickness_m(geometry, 100.0, 1.5)
    slow = ehl_film_thickness_m(geometry, 100.0, 0.5)
    heavy = ehl_film_thickness_m(geometry, 500.0, 1.5)
    assert fast > slow > 0.0
    assert heavy < fast
    assert ehl_film_thickness_m(geometry, 0.0, 1.0) == float("inf")


def test_self_consistent_invalid_inputs_raise(geometry):
    with pytest.raises(ValueError):
        integrate_sealing_contact(geometry, revolutions=0)
    with pytest.raises(ValueError):
        integrate_sealing_contact(geometry, rotor_mass_kg=-1.0)
