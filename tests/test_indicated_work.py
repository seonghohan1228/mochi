"""Indicated work and power from the port-timed pressure-volume loop."""

import pytest

from mochi.chambers import CycleTrace, build_cycle_trace
from mochi.indicated_work import indicated_work_j
from mochi.kinematics import RotaryGeometry

# A coarse cached trace keeps the suite fast; the value is grid-insensitive to a
# few percent, so the assertions below use loose magnitude bounds plus exact
# structural relations.
TRACE_SAMPLES = 72
TRACE_STEP_M = 2.0e-4


@pytest.fixture(scope="module")
def trace() -> CycleTrace:
    return build_cycle_trace(
        RotaryGeometry.default(), samples=TRACE_SAMPLES, grid_step_m=TRACE_STEP_M
    )


def test_power_is_work_times_frequency(trace: CycleTrace) -> None:
    geometry = RotaryGeometry.default()

    work = indicated_work_j(geometry, trace=trace)

    assert work.frequency_hz == pytest.approx(geometry.frequency_hz)
    assert work.power_w == pytest.approx(work.net_work_j * geometry.frequency_hz)
    assert work.net_work_j == pytest.approx(work.compression_work_j + work.suction_work_j)


def test_compressor_work_input_is_positive_and_of_the_right_scale(trace: CycleTrace) -> None:
    """The gas receives net work, dominated by the compression chamber.

    R410A 0.82 -> 3.24 MPa at n = 1.07 over a ~22 cm3 swept volume gives an
    indicated power of order several hundred watts; the exact figure tracks the
    grid but the scale and sign are robust.
    """

    work = indicated_work_j(RotaryGeometry.default(), trace=trace)

    assert work.net_work_j > 0.0  # net work input to the gas
    assert work.compression_work_j > 0.0  # compression chamber absorbs work
    assert work.suction_work_j < 0.0  # the filling suction chamber returns some
    assert 15.0 < work.net_work_j < 35.0
    assert 450.0 < work.power_w < 1050.0


def test_indicated_work_converges_with_the_sample_count(trace: CycleTrace) -> None:
    geometry = RotaryGeometry.default()

    coarse = indicated_work_j(geometry, trace=trace, samples=180)
    fine = indicated_work_j(geometry, trace=trace, samples=720)

    assert coarse.power_w == pytest.approx(fine.power_w, rel=0.01)


def test_too_few_samples_are_rejected(trace: CycleTrace) -> None:
    with pytest.raises(ValueError, match="at least eight"):
        indicated_work_j(RotaryGeometry.default(), trace=trace, samples=4)
