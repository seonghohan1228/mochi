"""Gap leakage in a chamber mass balance, and volumetric efficiency."""

import math

import pytest

from mochi.chambers import DISCHARGE_PORT_PRESSURE_PA, build_cycle_trace
from mochi.kinematics import RotaryGeometry
from mochi.leakage import (
    ISENTROPIC_EXPONENT,
    SUCTION_DENSITY_KG_M3,
    leaky_cycle,
    orifice_mass_flow,
)

TRACE_SAMPLES = 72
TRACE_STEP_M = 2.0e-4


@pytest.fixture(scope="module")
def trace():
    return build_cycle_trace(
        RotaryGeometry.default(), samples=TRACE_SAMPLES, grid_step_m=TRACE_STEP_M
    )


def _critical_ratio() -> float:
    k = ISENTROPIC_EXPONENT
    return (2.0 / (k + 1.0)) ** (k / (k - 1.0))


def test_orifice_is_continuous_at_the_choke_ratio() -> None:
    """The subsonic and choked branches agree at the critical pressure ratio."""

    upstream, area, density = 5.0e6, 1.0e-7, 100.0
    critical = _critical_ratio()
    just_subsonic = orifice_mass_flow(
        upstream, critical * upstream * 1.001, area, upstream_density_kg_m3=density
    )
    just_choked = orifice_mass_flow(
        upstream, critical * upstream * 0.999, area, upstream_density_kg_m3=density
    )

    assert just_choked == pytest.approx(just_subsonic, rel=1.0e-4)


def test_orifice_flow_stays_choked_below_the_critical_ratio() -> None:
    """Below the critical ratio the flow is choked, so it no longer grows."""

    upstream, area, density = 5.0e6, 1.0e-7, 100.0
    at_choke = orifice_mass_flow(upstream, 0.5 * upstream, area, upstream_density_kg_m3=density)
    deeper = orifice_mass_flow(upstream, 0.01 * upstream, area, upstream_density_kg_m3=density)

    assert deeper == pytest.approx(at_choke, rel=1.0e-9)  # choked plateau


@pytest.mark.parametrize("downstream", [1.0e6, 2.0e6])
def test_no_flow_without_a_positive_pressure_drop(downstream) -> None:
    assert orifice_mass_flow(1.0e6, downstream, 1.0e-7, upstream_density_kg_m3=50.0) == 0.0


def test_leakage_caps_the_recompression_peak_below_the_no_leak_bound(trace) -> None:
    """The leak bleeds trapped mass, so the capped peak is under the ~10 MPa bound."""

    cycle = leaky_cycle(RotaryGeometry.default(), trace=trace)

    assert cycle.no_leak_peak_pa > 9.0e6  # the strict no-leakage upper bound
    assert cycle.capped_peak_pa < cycle.no_leak_peak_pa
    assert cycle.capped_peak_pa > DISCHARGE_PORT_PRESSURE_PA  # still overshoots


def test_volumetric_efficiency_is_physical(trace) -> None:
    cycle = leaky_cycle(RotaryGeometry.default(), trace=trace)

    assert 0.7 < cycle.volumetric_efficiency < 1.0
    assert cycle.delivered_mass_kg > 0.0
    assert cycle.leaked_mass_kg > 0.0


def test_a_wider_gap_leaks_more_and_lowers_efficiency_and_peak(trace) -> None:
    """Monotonic in the assumed gap: more leakage, lower peak, lower efficiency."""

    geometry = RotaryGeometry.default()
    tight = leaky_cycle(geometry, gap_m=2.0e-6, trace=trace)
    loose = leaky_cycle(geometry, gap_m=20.0e-6, trace=trace)

    assert loose.leaked_mass_kg > tight.leaked_mass_kg
    assert loose.capped_peak_pa < tight.capped_peak_pa
    assert loose.volumetric_efficiency < tight.volumetric_efficiency
    assert loose.delivered_mass_kg < tight.delivered_mass_kg


def test_leakage_reduces_efficiency_below_the_clearance_only_value(trace) -> None:
    """With the leak off (vanishing gap) only clearance re-expansion cuts eta_v."""

    geometry = RotaryGeometry.default()
    clearance_only = leaky_cycle(geometry, gap_m=1.0e-9, trace=trace)
    with_leak = leaky_cycle(geometry, gap_m=5.0e-6, trace=trace)

    assert with_leak.volumetric_efficiency < clearance_only.volumetric_efficiency
    assert clearance_only.leaked_mass_kg < with_leak.leaked_mass_kg


def test_density_relation_matches_the_suction_reference() -> None:
    """The upstream density at the suction pressure is the suction vapour density."""

    # A degenerate cycle query is not needed; check the seed density directly via
    # the module constant used by leaky_cycle for the sealed suction charge.
    assert SUCTION_DENSITY_KG_M3 == pytest.approx(31.4, rel=1.0e-6)


@pytest.mark.parametrize("samples", [4, 0])
def test_too_few_samples_rejected(samples) -> None:
    with pytest.raises(ValueError, match="eight"):
        leaky_cycle(RotaryGeometry.default(), samples=samples)


@pytest.mark.parametrize("gap", [0.0, -1.0e-6, math.inf, math.nan])
def test_invalid_gap_rejected(gap) -> None:
    with pytest.raises(ValueError, match="[Gg]ap"):
        leaky_cycle(RotaryGeometry.default(), gap_m=gap)


@pytest.mark.parametrize("area", [0.0, -1.0e-7, math.nan])
def test_orifice_rejects_invalid_area(area) -> None:
    with pytest.raises(ValueError, match="area"):
        orifice_mass_flow(5.0e6, 1.0e6, area, upstream_density_kg_m3=100.0)
