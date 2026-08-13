"""Infinitely-long (Sommerfeld) journal bearing and the arc-film cross-check."""

import math

import numpy as np
import pytest

from mochi.arc_film import arc_film_force
from mochi.long_bearing import (
    half_sommerfeld_attitude_rad,
    long_bearing_load,
    sommerfeld_pressure,
)

# A representative long bearing (bush-scale radius, L/D > 1) at a plain entrainment.
RADIUS_M = 7.97e-3
LENGTH_M = 20.0e-3
CLEARANCE_M = 30.0e-6
ENTRAINMENT = 50.0
VISCOSITY = 0.010
EPS_SWEEP = (0.2, 0.4, 0.6, 0.8, 0.9)


# ---------------------------------------------------------------------------
# The analytical model is internally consistent
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("eps", EPS_SWEEP)
def test_sommerfeld_pressure_satisfies_the_reynolds_equation(eps):
    """d/dtheta(h^3 dp/dtheta) = 6 mu U R dh/dtheta for the closed-form pressure."""

    theta = np.linspace(0.05, math.pi - 0.05, 20001)
    d = theta[1] - theta[0]
    p = sommerfeld_pressure(
        eps,
        theta,
        radius_m=RADIUS_M,
        clearance_m=CLEARANCE_M,
        entrainment_speed_rad_s=ENTRAINMENT,
        viscosity_pa_s=VISCOSITY,
    )
    h = CLEARANCE_M * (1.0 + eps * np.cos(theta))
    flux = h**3 * np.gradient(p, d)
    lhs = np.gradient(flux, d)
    surface_speed = 2.0 * RADIUS_M * ENTRAINMENT
    rhs = 6.0 * VISCOSITY * surface_speed * RADIUS_M * np.gradient(h, d)
    # Compare on the interior (away from the differencing edges), relative to the scale.
    interior = slice(50, -50)
    scale = np.max(np.abs(rhs[interior]))
    assert np.max(np.abs(lhs[interior] - rhs[interior])) < 1e-3 * scale


@pytest.mark.parametrize("eps", EPS_SWEEP)
def test_full_sommerfeld_radial_load_vanishes(eps):
    """Fully flooded, the load is purely tangential (attitude 90 deg)."""

    load = long_bearing_load(
        eps,
        ENTRAINMENT,
        radius_m=RADIUS_M,
        length_m=LENGTH_M,
        clearance_m=CLEARANCE_M,
        viscosity_pa_s=VISCOSITY,
        condition="full",
    )
    assert abs(load.radial_force_n) < 1e-6 * abs(load.tangential_force_n)
    assert load.attitude_angle_rad == pytest.approx(math.pi / 2, abs=1e-6)


@pytest.mark.parametrize("eps", EPS_SWEEP)
def test_half_sommerfeld_attitude_matches_the_invariant(eps):
    """The Gumbel attitude equals atan(pi sqrt(1-eps^2)/(2 eps)), independent of scale."""

    load = long_bearing_load(
        eps,
        ENTRAINMENT,
        radius_m=RADIUS_M,
        length_m=LENGTH_M,
        clearance_m=CLEARANCE_M,
        viscosity_pa_s=VISCOSITY,
        condition="half",
    )
    assert load.attitude_angle_rad == pytest.approx(half_sommerfeld_attitude_rad(eps), rel=1e-4)


# ---------------------------------------------------------------------------
# arc_film reproduces the analytical long bearing (the validation)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("eps", EPS_SWEEP)
def test_arc_film_matches_long_bearing_sommerfeld(eps):
    """The bush curved-film solver equals the closed-form long bearing on its arc.

    The exact Sommerfeld pressure is zero at the maximum-film point (theta = 0) and at
    the minimum-film point (theta = pi), so over the converging half 0..pi the film is
    a Dirichlet problem with p = 0 ends -- exactly what ``arc_film`` solves -- and the
    pressure stays positive (Gumbel cavitation inactive). ``arc_film`` must therefore
    reproduce the closed form to its discretisation error. Film ``h = c(1 + eps cos
    beta)`` needs the piece eccentricity ``e_x = -eps c`` (min film at beta = pi).
    """

    arc = arc_film_force(
        -eps * CLEARANCE_M,
        0.0,
        0.0,
        0.0,
        ENTRAINMENT,
        arc_center_rad=math.pi / 2,
        arc_half_span_rad=math.pi / 2,  # spans [0, pi], the converging half
        radius_m=RADIUS_M,
        length_m=LENGTH_M,
        clearance_m=CLEARANCE_M,
        viscosity_pa_s=VISCOSITY,
        n_beta=4001,
    )
    reference = long_bearing_load(
        eps,
        ENTRAINMENT,
        radius_m=RADIUS_M,
        length_m=LENGTH_M,
        clearance_m=CLEARANCE_M,
        viscosity_pa_s=VISCOSITY,
        condition="half",
    )
    assert arc.force_x_n == pytest.approx(reference.radial_force_n, rel=1e-4)
    assert arc.force_y_n == pytest.approx(reference.tangential_force_n, rel=1e-4)
    assert arc.max_pressure_pa == pytest.approx(reference.max_pressure_pa, rel=1e-3)
    assert arc.min_film_thickness_m == pytest.approx(reference.min_film_thickness_m, rel=1e-6)


def test_arc_film_matches_long_bearing_when_rotated():
    """The agreement is orientation-independent (rotate the eccentricity by an angle)."""

    eps = 0.6
    phi = 0.7  # rotate the max-film axis to beta = phi
    arc = arc_film_force(
        -eps * CLEARANCE_M * math.cos(phi),
        -eps * CLEARANCE_M * math.sin(phi),
        0.0,
        0.0,
        ENTRAINMENT,
        arc_center_rad=phi + math.pi / 2,  # converging half runs from phi to phi + pi
        arc_half_span_rad=math.pi / 2,
        radius_m=RADIUS_M,
        length_m=LENGTH_M,
        clearance_m=CLEARANCE_M,
        viscosity_pa_s=VISCOSITY,
        n_beta=4001,
    )
    reference = long_bearing_load(
        eps,
        ENTRAINMENT,
        radius_m=RADIUS_M,
        length_m=LENGTH_M,
        clearance_m=CLEARANCE_M,
        viscosity_pa_s=VISCOSITY,
        condition="half",
    )
    # Compare magnitudes (the resultant is frame-independent).
    assert math.hypot(arc.force_x_n, arc.force_y_n) == pytest.approx(
        reference.magnitude_n, rel=1e-4
    )


def test_invalid_inputs_raise():
    with pytest.raises(ValueError):
        long_bearing_load(
            1.0, ENTRAINMENT, radius_m=RADIUS_M, length_m=LENGTH_M, clearance_m=CLEARANCE_M
        )
    with pytest.raises(ValueError):
        long_bearing_load(
            0.5,
            ENTRAINMENT,
            radius_m=RADIUS_M,
            length_m=LENGTH_M,
            clearance_m=CLEARANCE_M,
            condition="bogus",
        )
