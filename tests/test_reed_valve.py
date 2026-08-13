"""Discharge reed valve on the chamber mass balance (quasi-static check valve)."""

import math

import pytest

from mochi.chambers import DISCHARGE_PORT_PRESSURE_PA, build_cycle_trace
from mochi.kinematics import RotaryGeometry
from mochi.reed_valve import valved_cycle

TRACE_SAMPLES = 72
TRACE_STEP_M = 2.0e-4


@pytest.fixture(scope="module")
def trace():
    return build_cycle_trace(
        RotaryGeometry.default(), samples=TRACE_SAMPLES, grid_step_m=TRACE_STEP_M
    )


def test_finite_port_area_makes_the_discharge_overpressure(trace) -> None:
    """The valve cannot vent instantly, so delivery runs above the line pressure."""

    cycle = valved_cycle(RotaryGeometry.default(), trace=trace)

    assert cycle.delivery_peak_pa > DISCHARGE_PORT_PRESSURE_PA  # overpressure
    assert cycle.overpressure_work_j > 0.0
    assert cycle.overpressure_power_w > 0.0


def test_a_wider_valve_lowers_the_overpressure_and_loss(trace) -> None:
    """More flow area lets the valve keep pace, so the overpressure shrinks."""

    geometry = RotaryGeometry.default()
    tight = valved_cycle(geometry, valve_area_scale=0.5, trace=trace)
    wide = valved_cycle(geometry, valve_area_scale=4.0, trace=trace)

    assert wide.delivery_peak_pa < tight.delivery_peak_pa
    assert wide.overpressure_power_w < tight.overpressure_power_w


def test_overpressure_floor_is_the_reed_opening_pressure(trace) -> None:
    """Even a very large valve cannot deliver below the reed opening pressure."""

    geometry = RotaryGeometry.default()
    huge = valved_cycle(geometry, valve_area_scale=50.0, trace=trace)

    opening = DISCHARGE_PORT_PRESSURE_PA * 1.05
    assert huge.delivery_peak_pa == pytest.approx(opening, rel=0.02)


def test_valve_indicated_power_exceeds_the_baseline(trace) -> None:
    """The overpressure adds P-V area, so the valve indicated power > the ideal."""

    cycle = valved_cycle(RotaryGeometry.default(), trace=trace)

    assert cycle.valve_indicated_work_j == pytest.approx(
        cycle.baseline_indicated_work_j + cycle.overpressure_work_j, rel=1.0e-12
    )
    assert cycle.valve_indicated_power_w > cycle.baseline_indicated_work_j * cycle.frequency_hz
    assert 770.0 < cycle.valve_indicated_power_w < 810.0  # ~789 W = 738 + 51


def test_the_valve_cannot_cap_the_recompression_spike(trace) -> None:
    """After the geometric port closes the area is zero, so the spike remains."""

    cycle = valved_cycle(RotaryGeometry.default(), trace=trace)

    assert cycle.recompression_peak_pa > 8.0e6  # still the leakage-limited spike


def test_volumetric_efficiency_and_masses_are_physical(trace) -> None:
    cycle = valved_cycle(RotaryGeometry.default(), trace=trace)

    assert 0.8 < cycle.volumetric_efficiency < 1.0
    assert cycle.delivered_mass_kg > 0.0
    assert cycle.leaked_mass_kg > 0.0


@pytest.mark.parametrize("samples", [4, 0])
def test_too_few_samples_rejected(samples) -> None:
    with pytest.raises(ValueError, match="eight"):
        valved_cycle(RotaryGeometry.default(), samples=samples)


@pytest.mark.parametrize("gap", [0.0, -1.0e-6, math.inf, math.nan])
def test_invalid_gap_rejected(gap) -> None:
    with pytest.raises(ValueError, match="[Gg]ap"):
        valved_cycle(RotaryGeometry.default(), gap_m=gap)
