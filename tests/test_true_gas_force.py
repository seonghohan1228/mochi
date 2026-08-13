"""Mouth-aware true-geometry gas force and shaft torque (PHYSICS.md 4.5)."""

import math

import pytest

from mochi.chambers import build_cycle_trace
from mochi.gas_force import gas_load
from mochi.indicated_work import indicated_work_j
from mochi.kinematics import RotaryGeometry
from mochi.true_gas_force import (
    peak_rotor_force_n,
    true_gas_load,
    true_gas_torque_work_j,
)

TRACE_SAMPLES = 72
TRACE_STEP_M = 2.0e-4
CROSS_CHECK_SAMPLES = 180


@pytest.fixture(scope="module")
def trace():
    return build_cycle_trace(
        RotaryGeometry.default(), samples=TRACE_SAMPLES, grid_step_m=TRACE_STEP_M
    )


def test_gas_moment_about_rotor_centre_is_nonzero(trace) -> None:
    """The mouth asymmetry offsets the centre of pressure -> a moment about O_b (Section 4.13)."""

    geometry = RotaryGeometry.default()
    load = true_gas_load(geometry, math.radians(226.0), trace=trace)
    # the true resultant does not pass through O_r, so a finite moment exists
    assert abs(load.rotor_moment_about_center_nm) > 1.0
    assert math.isfinite(load.rotor_moment_about_center_nm)


def test_torque_work_equals_the_indicated_work(trace) -> None:
    """The true shaft-torque integral is the true indicated work, above circular."""

    geometry = RotaryGeometry.default()
    work = true_gas_torque_work_j(geometry, samples=CROSS_CHECK_SAMPLES, trace=trace)
    indicated = indicated_work_j(geometry, samples=CROSS_CHECK_SAMPLES, trace=trace)

    # oint T dtheta = -oint p dV = W, evaluated seam-safely (Section 4.5).
    assert work.torque_work_j == pytest.approx(indicated.net_work_j, rel=1.0e-9)
    assert work.power_w == pytest.approx(indicated.power_w, rel=1.0e-9)
    # The rotor mouth lifts the circular-rotor closed form (~3%, 715 -> 738 W).
    assert work.power_w > work.circular_power_w
    assert work.power_w / work.circular_power_w == pytest.approx(1.03, abs=0.02)


@pytest.mark.parametrize("theta_deg", [200.0, 235.0, 300.0])
def test_true_force_is_below_the_circular_closed_form(theta_deg, trace) -> None:
    """The mouth cavity load lowers the net rotor force below the OD-disc form."""

    geometry = RotaryGeometry.default()
    theta = math.radians(theta_deg)
    true_load = true_gas_load(geometry, theta, trace=trace)
    circular = gas_load(geometry, theta, trace=trace)

    circular_mag = math.hypot(*circular.rotor_force_n)
    assert true_load.rotor_force_mag_n < circular_mag
    assert true_load.rotor_force_mag_n == pytest.approx(
        math.hypot(*true_load.rotor_force_n), rel=1.0e-12
    )


def test_peak_rotor_force_is_physical_and_below_circular(trace) -> None:
    """The peak reaction is kN-scale but ~20% under the circular peak (~3.3 kN)."""

    geometry = RotaryGeometry.default()
    peak = peak_rotor_force_n(geometry, trace=trace)

    assert 1.0e3 < peak < 3.0e3
    # Well below the circular closed-form peak of Section 4.6.
    assert peak < 3.0e3


def test_seal_over_window_is_rejected() -> None:
    with pytest.raises(ValueError, match="seal-over"):
        true_gas_load(RotaryGeometry.default(), 0.0)


@pytest.mark.parametrize("angle", [math.inf, math.nan])
def test_non_finite_angle_is_rejected(angle) -> None:
    with pytest.raises(ValueError, match="finite"):
        true_gas_load(RotaryGeometry.default(), angle)


def test_torque_work_rejects_too_few_samples() -> None:
    with pytest.raises(ValueError, match="eight"):
        true_gas_torque_work_j(RotaryGeometry.default(), samples=4)


def test_peak_force_rejects_too_few_samples() -> None:
    with pytest.raises(ValueError, match="eight"):
        peak_rotor_force_n(RotaryGeometry.default(), samples=4)
