"""Crank-pin bearing reaction and drive torque (journal-lumped, bush friction)."""

import math

import pytest

from mochi.bearing_load import mechanism_load, shaft_work_j
from mochi.chambers import build_cycle_trace
from mochi.gas_force import gas_load, gas_torque_work_j
from mochi.kinematics import RotaryGeometry

TRACE_SAMPLES = 72
TRACE_STEP_M = 2.0e-4
CROSS_CHECK_SAMPLES = 180


@pytest.fixture(scope="module")
def trace():
    return build_cycle_trace(
        RotaryGeometry.default(), samples=TRACE_SAMPLES, grid_step_m=TRACE_STEP_M
    )


@pytest.fixture(scope="module")
def work(trace):
    return shaft_work_j(RotaryGeometry.default(), samples=CROSS_CHECK_SAMPLES, trace=trace)


@pytest.mark.parametrize("theta_deg", [45.0, 120.0, 230.0, 300.0])
def test_journal_reaction_is_minus_the_gas_force(theta_deg, trace) -> None:
    """The journal carries the whole gas load: R_j = -F_gas (Section 4.6)."""

    geometry = RotaryGeometry.default()
    theta = math.radians(theta_deg)
    load = mechanism_load(geometry, theta, trace=trace)
    gas = gas_load(geometry, theta, trace=trace)

    assert load.journal_force_n[0] == pytest.approx(-gas.rotor_force_n[0], rel=1.0e-12)
    assert load.journal_force_n[1] == pytest.approx(-gas.rotor_force_n[1], rel=1.0e-12)
    assert load.journal_force_mag_n == pytest.approx(math.hypot(*gas.rotor_force_n), rel=1.0e-12)


def test_drive_torque_splits_into_gas_plus_friction(trace) -> None:
    """T_drive = T_gas + (bush + journal) friction, all positive dissipation."""

    load = mechanism_load(RotaryGeometry.default(), math.radians(230.0), trace=trace)

    assert load.friction_torque_nm == pytest.approx(
        load.bush_friction_torque_nm + load.journal_friction_torque_nm, rel=1.0e-12
    )
    assert load.drive_torque_nm == pytest.approx(
        load.gas_torque_nm + load.friction_torque_nm, rel=1.0e-12
    )
    assert load.bush_friction_torque_nm > 0.0
    assert load.journal_friction_torque_nm > 0.0


def test_drive_work_splits_exactly_into_gas_and_friction(work) -> None:
    """The integrated drive work is exactly gas work plus friction work."""

    assert work.drive_work_j == pytest.approx(work.gas_work_j + work.friction_work_j, rel=1.0e-12)
    assert work.friction_work_j > 0.0


def test_gas_work_matches_the_gas_torque_rung(work, trace) -> None:
    """The gas-torque part reproduces the Section 4.5 gas-torque work."""

    gas = gas_torque_work_j(RotaryGeometry.default(), samples=CROSS_CHECK_SAMPLES, trace=trace)

    assert work.gas_work_j == pytest.approx(gas.torque_work_j, rel=1.0e-9)


def test_shaft_power_is_indicated_plus_friction(work) -> None:
    """Shaft power = baseline indicated + friction (reed valve not propagated)."""

    # Exact by construction (Section 3.5 baseline indicated + friction, Section
    # 4.6/4.7); the reed-valve overpressure (Section 3.8) is a separate term.
    assert work.shaft_power_w == pytest.approx(
        work.indicated_power_w + work.friction_power_w, rel=1.0e-12
    )
    # The shaft supplies the indicated gas work plus the losses, above indicated.
    assert work.shaft_power_w > work.indicated_power_w
    # Friction splits into bush (~0.2 W) and journal (~9 W), journal dominant.
    assert work.friction_power_w == pytest.approx(
        work.bush_friction_power_w + work.journal_friction_power_w, rel=1.0e-12
    )
    assert 0.15 < work.bush_friction_power_w < 0.30
    assert 8.0 < work.journal_friction_power_w < 20.0
    assert work.journal_friction_power_w > 10.0 * work.bush_friction_power_w


def test_peak_journal_load_is_positive_and_physical(work) -> None:
    """A single-cylinder R410A machine of this size carries a kN-scale bearing load."""

    assert 1.0e3 < work.peak_journal_force_n < 1.0e4


def test_seal_over_window_is_rejected() -> None:
    with pytest.raises(ValueError, match="seal-over"):
        mechanism_load(RotaryGeometry.default(), 0.0)


@pytest.mark.parametrize("angle", [math.inf, math.nan])
def test_non_finite_angle_is_rejected(angle) -> None:
    with pytest.raises(ValueError, match="finite"):
        mechanism_load(RotaryGeometry.default(), angle)


def test_shaft_work_rejects_too_few_samples() -> None:
    with pytest.raises(ValueError, match="eight"):
        shaft_work_j(RotaryGeometry.default(), samples=4)
