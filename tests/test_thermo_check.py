"""Thermodynamic cross-check of the indicated work (CoolProp isentropic route)."""

import pytest

from mochi.chambers import build_cycle_trace
from mochi.kinematics import RotaryGeometry
from mochi.thermo_check import isentropic_cross_check

TRACE_SAMPLES = 72
TRACE_STEP_M = 2.0e-4


@pytest.fixture(scope="module")
def trace():
    return build_cycle_trace(
        RotaryGeometry.default(), samples=TRACE_SAMPLES, grid_step_m=TRACE_STEP_M
    )


def test_indicated_work_matches_the_isentropic_route(trace) -> None:
    """Route A (P-V loop) and Route B (m * dh_s) agree within a couple percent."""

    check = isentropic_cross_check(RotaryGeometry.default(), trace=trace)

    assert check.relative_error < 0.02
    assert check.indicated_power_w == pytest.approx(check.isentropic_power_w, rel=0.02)


def test_specific_work_is_the_isentropic_enthalpy_rise(trace) -> None:
    """W/m equals the CoolProp isentropic enthalpy rise for the near-isentropic cycle."""

    check = isentropic_cross_check(RotaryGeometry.default(), trace=trace)

    assert check.specific_work_j_kg == pytest.approx(check.isentropic_enthalpy_rise_j_kg, rel=0.02)
    assert 30.0e3 < check.specific_work_j_kg < 45.0e3
    assert check.delivered_mass_kg > 0.0
    assert check.indicated_work_j == pytest.approx(
        check.specific_work_j_kg * check.delivered_mass_kg, rel=1.0e-9
    )


def test_non_positive_enthalpy_rise_is_rejected() -> None:
    with pytest.raises(ValueError, match="enthalpy"):
        isentropic_cross_check(RotaryGeometry.default(), enthalpy_rise_j_kg=0.0)
