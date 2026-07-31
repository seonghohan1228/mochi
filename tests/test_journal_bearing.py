"""Crank-pin journal bearing friction (Petroff concentric hydrodynamic model)."""

import math

import pytest

from mochi.journal_bearing import (
    JOURNAL_CLEARANCE_M,
    JOURNAL_LENGTH_M,
    JOURNAL_RADIUS_M,
    LUBRICANT_VISCOSITY_PA_S,
    journal_relative_speed_rad_s,
    petroff_friction,
)
from mochi.kinematics import RotaryGeometry


def _petroff_torque(speed_rad_s: float) -> float:
    return (
        2.0
        * math.pi
        * LUBRICANT_VISCOSITY_PA_S
        * speed_rad_s
        * JOURNAL_RADIUS_M**3
        * JOURNAL_LENGTH_M
        / JOURNAL_CLEARANCE_M
    )


@pytest.mark.parametrize("theta_deg", [10.0, 90.0, 180.0, 270.0, 340.0])
def test_torque_matches_the_petroff_formula(theta_deg) -> None:
    """T_j = 2 pi mu |omega_rel| r^3 L / c, evaluated at the model geometry."""

    geometry = RotaryGeometry.default()
    theta = math.radians(theta_deg)
    friction = petroff_friction(geometry, theta)
    speed = abs(journal_relative_speed_rad_s(geometry, theta))

    assert friction.friction_torque_nm == pytest.approx(_petroff_torque(speed), rel=1.0e-12)


def test_friction_power_is_torque_times_relative_speed() -> None:
    geometry = RotaryGeometry.default()
    theta = math.radians(180.0)
    friction = petroff_friction(geometry, theta)

    assert friction.friction_power_w == pytest.approx(
        friction.friction_torque_nm * abs(friction.relative_speed_rad_s), rel=1.0e-12
    )


def test_relative_speed_is_the_shaft_speed_minus_the_rotor_swing() -> None:
    """omega_rel = omega*(1 - dphi/dtheta); it stays near, and below 2x, the shaft speed."""

    geometry = RotaryGeometry.default()
    shaft = geometry.angular_speed_rad_s

    speeds = [journal_relative_speed_rad_s(geometry, math.radians(d)) for d in range(0, 360, 5)]
    # The rotor only swings +/-10.37 deg, so the relative speed brackets the shaft
    # speed and never reverses (150-225 rad/s at the default 188.5 rad/s).
    assert min(speeds) > 150.0
    assert max(speeds) < 225.0
    # At the swing extremes (theta = 90, 270 deg) the rotor spin is zero, so the
    # relative speed equals the shaft speed.
    assert journal_relative_speed_rad_s(geometry, math.radians(90.0)) == pytest.approx(
        shaft, rel=1.0e-3
    )


def test_torque_scales_with_viscosity_and_inverse_clearance() -> None:
    geometry = RotaryGeometry.default()
    theta = math.radians(120.0)
    base = petroff_friction(geometry, theta).friction_torque_nm

    twice_viscosity = petroff_friction(
        geometry, theta, viscosity_pa_s=2.0 * LUBRICANT_VISCOSITY_PA_S
    )
    half_clearance = petroff_friction(geometry, theta, clearance_m=0.5 * JOURNAL_CLEARANCE_M)

    assert twice_viscosity.friction_torque_nm == pytest.approx(2.0 * base, rel=1.0e-12)
    assert half_clearance.friction_torque_nm == pytest.approx(2.0 * base, rel=1.0e-12)


def test_sommerfeld_number_is_positive_under_load_and_infinite_without() -> None:
    geometry = RotaryGeometry.default()
    theta = math.radians(200.0)

    loaded = petroff_friction(geometry, theta, load_n=3000.0)
    unloaded = petroff_friction(geometry, theta)

    assert 0.0 < loaded.sommerfeld_number < 1.0  # heavily loaded compressor bearing
    assert math.isinf(unloaded.sommerfeld_number)


def test_journal_power_is_a_few_percent_of_indicated() -> None:
    """The assumed geometry gives an order-10 W loss, dwarfing the bush film."""

    geometry = RotaryGeometry.default()
    powers = [
        petroff_friction(geometry, math.radians(d)).friction_power_w for d in range(0, 360, 5)
    ]
    mean_power = sum(powers) / len(powers)

    assert 8.0 < mean_power < 20.0


@pytest.mark.parametrize("angle", [math.inf, math.nan])
def test_non_finite_angle_is_rejected(angle) -> None:
    with pytest.raises(ValueError, match="finite"):
        petroff_friction(RotaryGeometry.default(), angle)


@pytest.mark.parametrize("clearance", [0.0, -1.0e-6, math.inf, math.nan])
def test_invalid_clearance_is_rejected(clearance) -> None:
    with pytest.raises(ValueError, match="clearance"):
        petroff_friction(RotaryGeometry.default(), 0.0, clearance_m=clearance)
