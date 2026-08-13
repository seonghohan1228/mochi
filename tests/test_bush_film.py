"""Swing-bush Couette + Poiseuille lubrication-film model."""

import math

import pytest

from mochi.bush_film import (
    AXIAL_CLEARANCE_M,
    LUBRICANT_VISCOSITY_PA_S,
    curved_slide_velocity,
    film_state,
    film_thicknesses_m,
    flat_contact_length_m,
    flat_slide_velocity,
    friction_power_cycle_w,
    full_flat_contact_length_m,
)
from mochi.chamber_volume import SwingBush
from mochi.chambers import DISCHARGE_PORT_PRESSURE_PA, SUCTION_PORT_PRESSURE_PA, build_cycle_trace
from mochi.kinematics import MM, RotaryGeometry

TRACE_SAMPLES = 72
TRACE_STEP_M = 2.0e-4


@pytest.fixture(scope="module")
def trace():
    return build_cycle_trace(
        RotaryGeometry.default(), samples=TRACE_SAMPLES, grid_step_m=TRACE_STEP_M
    )


def test_flat_contact_shortens_near_bdc_from_the_vane_tip_round() -> None:
    """The R1.5 vane-tip round trims the bush flat contact near BDC, full elsewhere."""

    geometry = RotaryGeometry.default()
    full = full_flat_contact_length_m(SwingBush())

    assert full == pytest.approx(11.94 * MM, abs=0.01 * MM)
    # Away from BDC the bush flat sits fully on the straight flank.
    assert flat_contact_length_m(geometry, math.radians(0.0)) == pytest.approx(full, rel=1.0e-9)
    assert flat_contact_length_m(geometry, math.radians(90.0)) == pytest.approx(full, rel=1.0e-9)
    # At BDC the flat overhangs the tip round, shortening to ~11.47 mm.
    bdc = flat_contact_length_m(geometry, math.radians(180.0))
    assert bdc < full
    assert bdc == pytest.approx(11.47 * MM, abs=0.03 * MM)
    # Monotone: nearer BDC is not longer.
    assert bdc <= flat_contact_length_m(geometry, math.radians(165.0))


def test_both_films_are_ten_microns() -> None:
    flat_film, curved_film = film_thicknesses_m(RotaryGeometry.default(), SwingBush())

    assert flat_film == pytest.approx(0.010 * MM, rel=1.0e-6)
    assert curved_film == pytest.approx(0.010 * MM, rel=1.0e-6)


def test_sliding_velocities_match_the_supplied_peak_demands() -> None:
    """Translation peaks near 0.86 m/s, the curved swing near 0.27 m/s."""

    geometry = RotaryGeometry.default()
    angles = [math.radians(d) for d in range(360)]

    flat_peak = max(abs(flat_slide_velocity(geometry, a)) for a in angles)
    curved_peak = max(abs(curved_slide_velocity(geometry, a)) for a in angles)

    assert flat_peak == pytest.approx(0.86, rel=0.03)
    assert curved_peak == pytest.approx(0.27, rel=0.05)
    # Translation peak is the eccentricity times shaft speed, near theta = 90 deg.
    assert abs(flat_slide_velocity(geometry, math.radians(90.0))) == pytest.approx(
        geometry.eccentricity_m * geometry.angular_speed_rad_s, rel=1.0e-3
    )


def test_shear_stress_is_mu_u_over_h() -> None:
    geometry = RotaryGeometry.default()
    theta = math.radians(90.0)
    velocity = flat_slide_velocity(geometry, theta)

    face = film_state(
        geometry, theta, in_chamber_pressure_pa=1.0e6, out_chamber_pressure_pa=1.0e6
    ).in_piece.flat

    assert face.shear_stress_pa == pytest.approx(
        LUBRICANT_VISCOSITY_PA_S * velocity / (0.010 * MM), rel=1.0e-9
    )
    assert face.friction_power_w >= 0.0


def test_film_pressure_is_linear_between_the_end_pressures() -> None:
    face = film_state(
        RotaryGeometry.default(),
        math.radians(200.0),
        in_chamber_pressure_pa=2.0e6,
        out_chamber_pressure_pa=2.0e6,
        recess_pressure_pa=4.0e6,
    ).in_piece.flat

    assert face.pressure_at(0.0) == pytest.approx(face.inlet_pressure_pa)
    assert face.pressure_at(1.0) == pytest.approx(face.outlet_pressure_pa)
    assert face.pressure_at(0.5) == pytest.approx(
        0.5 * (face.inlet_pressure_pa + face.outlet_pressure_pa)
    )
    assert face.inlet_pressure_pa == pytest.approx(2.0e6)
    assert face.outlet_pressure_pa == pytest.approx(4.0e6)


def test_leakage_reduces_to_pure_couette_when_the_ends_are_equal() -> None:
    """With no pressure difference the Poiseuille term vanishes: q = U h/2 · H."""

    geometry = RotaryGeometry.default()
    theta = math.radians(90.0)
    height = geometry.cylinder_height_m - 2.0 * AXIAL_CLEARANCE_M

    face = film_state(
        geometry,
        theta,
        in_chamber_pressure_pa=1.0e6,
        out_chamber_pressure_pa=1.0e6,
        recess_pressure_pa=1.0e6,
    ).in_piece.flat

    assert face.leakage_m3_s == pytest.approx(
        face.slide_velocity_m_s * (0.010 * MM) / 2.0 * height, rel=1.0e-6
    )


def test_default_boundary_conditions_use_the_chamber_pressures(trace) -> None:
    geometry = RotaryGeometry.default()

    state = film_state(geometry, math.radians(270.0), trace=trace)

    # IN piece sees the suction chamber, OUT piece the compression chamber.
    assert state.in_piece.flat.inlet_pressure_pa == pytest.approx(SUCTION_PORT_PRESSURE_PA)
    assert state.out_piece.flat.inlet_pressure_pa == pytest.approx(
        DISCHARGE_PORT_PRESSURE_PA, rel=1.0e-3
    )  # at 270 deg the compression chamber is delivering at the discharge pressure
    assert state.in_piece.curved.outlet_pressure_pa == pytest.approx(4.0e6)


def test_cycle_mean_friction_power_is_positive_and_small(trace) -> None:
    power = friction_power_cycle_w(RotaryGeometry.default(), samples=72, trace=trace)

    assert 0.0 < power < 1.0  # order 0.1-0.3 W, tiny beside the indicated power


@pytest.mark.parametrize("viscosity", [0.0, -1.0, math.inf, math.nan])
def test_invalid_viscosity_is_rejected(viscosity: float) -> None:
    with pytest.raises(ValueError, match="[Vv]iscosity"):
        film_state(
            RotaryGeometry.default(),
            0.0,
            viscosity_pa_s=viscosity,
            in_chamber_pressure_pa=1e6,
            out_chamber_pressure_pa=1e6,
        )


@pytest.mark.parametrize("angle", [math.inf, math.nan])
def test_non_finite_angle_is_rejected(angle: float) -> None:
    with pytest.raises(ValueError, match="finite"):
        film_state(
            RotaryGeometry.default(), angle, in_chamber_pressure_pa=1e6, out_chamber_pressure_pa=1e6
        )
