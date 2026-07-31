"""Gas pressure force on the rotor and the net gas torque about the crank axis."""

import math

import pytest

from mochi.chambers import build_cycle_trace
from mochi.gas_force import (
    _arc_force,
    _foot_normal_angles,
    gas_load,
    gas_torque_work_j,
)
from mochi.indicated_work import indicated_work_j
from mochi.kinematics import RotaryGeometry

TRACE_SAMPLES = 72
TRACE_STEP_M = 2.0e-4
CROSS_CHECK_SAMPLES = 360


@pytest.fixture(scope="module")
def trace():
    return build_cycle_trace(
        RotaryGeometry.default(), samples=TRACE_SAMPLES, grid_step_m=TRACE_STEP_M
    )


def _boundary_integral(
    geometry: RotaryGeometry,
    crank_angle_rad: float,
    suction_pa: float,
    compression_pa: float,
    segments: int = 4000,
) -> tuple[tuple[float, float], float]:
    """Independent midpoint integral of ``p (-n) dl`` over the wetted rotor arc."""

    eccentricity = geometry.eccentricity_m
    rotor_radius = geometry.rotor_radius_m
    height = geometry.cylinder_height_m
    center_x = eccentricity * math.sin(crank_angle_rad)
    center_y = eccentricity * math.cos(crank_angle_rad)
    inlet_foot, outlet_foot = _foot_normal_angles(geometry, crank_angle_rad)

    force_x = force_y = torque = 0.0
    arcs = (
        (inlet_foot, crank_angle_rad, suction_pa),
        (crank_angle_rad, outlet_foot + 2.0 * math.pi, compression_pa),
    )
    for phi_start, phi_end, pressure in arcs:
        step = (phi_end - phi_start) / segments
        coefficient = -pressure * rotor_radius * height * step
        for index in range(segments):
            phi = phi_start + (index + 0.5) * step
            normal_x, normal_y = math.sin(phi), math.cos(phi)
            d_fx = coefficient * normal_x
            d_fy = coefficient * normal_y
            point_x = center_x + rotor_radius * normal_x
            point_y = center_y + rotor_radius * normal_y
            force_x += d_fx
            force_y += d_fy
            torque += point_x * d_fy - point_y * d_fx
    return (force_x, force_y), torque


def test_torque_work_matches_the_reference_pv_work(trace) -> None:
    """The keystone identity: oint T_gas d(theta) = -oint p dV on the same basis."""

    result = gas_torque_work_j(RotaryGeometry.default(), samples=CROSS_CHECK_SAMPLES, trace=trace)

    assert result.torque_work_j == pytest.approx(result.reference_pv_work_j, rel=1.0e-3)
    assert result.torque_work_j > 0.0  # compressor: net work into the gas


def test_power_is_near_the_indicated_power(trace) -> None:
    """Torque-integral power matches the Section 3.5 headline within the basis gap."""

    geometry = RotaryGeometry.default()
    torque = gas_torque_work_j(geometry, samples=CROSS_CHECK_SAMPLES, trace=trace)
    indicated = indicated_work_j(geometry, trace=trace)

    # The residual is the circular-rotor vs true-volume basis, not an error.
    assert torque.power_w == pytest.approx(indicated.power_w, rel=0.06)


@pytest.mark.parametrize("theta_deg", [45.0, 120.0, 200.0, 300.0])
def test_torque_equals_the_cross_product_of_center_and_force(theta_deg, trace) -> None:
    """T_gas = O_r x F_rotor exactly: the resultant passes through the rotor center."""

    geometry = RotaryGeometry.default()
    theta = math.radians(theta_deg)
    load = gas_load(geometry, theta, trace=trace)

    center_x = geometry.eccentricity_m * math.sin(theta)
    center_y = geometry.eccentricity_m * math.cos(theta)
    cross = center_x * load.rotor_force_n[1] - center_y * load.rotor_force_n[0]

    assert load.rotor_torque_nm == pytest.approx(cross, rel=1.0e-12, abs=1.0e-12)


@pytest.mark.parametrize("theta_deg", [45.0, 120.0, 200.0, 300.0])
def test_closed_form_force_matches_the_boundary_integral(theta_deg, trace) -> None:
    """The endpoint closed form agrees with an independent discretized integral."""

    geometry = RotaryGeometry.default()
    theta = math.radians(theta_deg)
    load = gas_load(geometry, theta, trace=trace)
    (force_x, force_y), torque = _boundary_integral(
        geometry, theta, load.suction_pressure_pa, load.compression_pressure_pa
    )

    assert load.rotor_force_n[0] == pytest.approx(force_x, rel=1.0e-3)
    assert load.rotor_force_n[1] == pytest.approx(force_y, rel=1.0e-3)
    assert load.rotor_torque_nm == pytest.approx(torque, rel=1.0e-3)


def test_uniform_pressure_on_the_full_circle_gives_zero_force() -> None:
    """A closed circle at uniform pressure carries no net force (the arc closure)."""

    geometry = RotaryGeometry.default()
    force_x, force_y = _arc_force(geometry, 2.5e6, 0.7, 0.7 + 2.0 * math.pi)

    assert force_x == pytest.approx(0.0, abs=1.0e-6)
    assert force_y == pytest.approx(0.0, abs=1.0e-6)


def test_discharge_pressure_pushes_the_rotor_toward_the_inlet(trace) -> None:
    """High compression pressure sits on the -x side, so the net force is +x."""

    geometry = RotaryGeometry.default()
    # Late compression / delivery: compression pressure far exceeds suction.
    load = gas_load(geometry, math.radians(300.0), trace=trace)

    assert load.compression_pressure_pa > load.suction_pressure_pa
    assert load.rotor_force_n[0] > 0.0


def test_seal_over_window_is_rejected() -> None:
    """No separated arcs exist while the contact hides under the vane."""

    with pytest.raises(ValueError, match="seal-over"):
        gas_load(RotaryGeometry.default(), 0.0)


@pytest.mark.parametrize("angle", [math.inf, math.nan])
def test_non_finite_angle_is_rejected(angle) -> None:
    with pytest.raises(ValueError, match="finite"):
        gas_load(RotaryGeometry.default(), angle)


def test_cross_check_rejects_too_few_samples() -> None:
    with pytest.raises(ValueError, match="eight"):
        gas_torque_work_j(RotaryGeometry.default(), samples=4)
