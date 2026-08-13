"""Nine-DOF rotor + two-piece swing-bush multibody dynamics (D2, Stage 5)."""

import math

import numpy as np
import pytest

from mochi.arc_film import arc_film_force
from mochi.bush_film import AXIAL_CLEARANCE_M, film_thicknesses_m, friction_power_cycle_w
from mochi.chamber_volume import SwingBush
from mochi.chambers import build_cycle_trace
from mochi.kinematics import RotaryGeometry
from mochi.mass_properties import bush_piece_mass_properties
from mochi.rotor_bush_dynamics import (
    BUSH_PIECE_MASS_KG,
    ROTOR_POLAR_INERTIA_KG_M2,
    _piece_inertia_about_centre,
    integrate_rotor_bush_orbit,
)
from mochi.slider_film import flat_slider_film

TRACE_SAMPLES = 72
TRACE_STEP_M = 2.0e-4


@pytest.fixture(scope="module")
def geometry():
    return RotaryGeometry.default()


@pytest.fixture(scope="module")
def trace(geometry):
    return build_cycle_trace(geometry, samples=TRACE_SAMPLES, grid_step_m=TRACE_STEP_M)


@pytest.fixture(scope="module")
def orbit(geometry, trace):
    # Coarse but converged enough for the reduction / boundedness checks; the fast
    # squeeze transient decays inside the first revolution. gas_film_boundary is off so
    # these check the hydrodynamic bush + seal dynamics in isolation (the gas boundary
    # condition, which drives the flat film toward contact, has its own test below).
    return integrate_rotor_bush_orbit(
        geometry,
        revolutions=2,
        samples=40,
        grid_samples=90,
        n_beta=61,
        n_s=41,
        trace=trace,
        gas_film_boundary=False,
    )


# ---------------------------------------------------------------------------
# Reference-film reduction (spec Section 8.2): zero deviation -> 10 um films
# ---------------------------------------------------------------------------


def test_reference_curved_film_reduces_to_ten_microns(geometry):
    """The reference piece eccentricity (piece_shift) gives the 10 um curved film.

    Concentric gap r_cut - r_b = 30 um; the 20 um shift toward the piece's own side
    makes the minimum curved film 10 um, matching :func:`film_thicknesses_m`.
    """

    bush = SwingBush()
    radius = bush.piece_outer_radius_m
    height = geometry.cylinder_height_m - 2.0 * AXIAL_CLEARANCE_M
    curved_gap = geometry.cutout_radius_m - radius
    _, expected_curved = film_thicknesses_m(geometry, bush)

    # IN piece: reference eccentricity is piece_shift along +x (its own side).
    force = arc_film_force(
        bush.piece_shift_m,
        0.0,
        0.0,
        0.0,
        20.0,  # some entrainment; the minimum-film geometry is independent of it
        arc_center_rad=0.0,
        arc_half_span_rad=bush.half_arc_rad(),
        radius_m=radius,
        length_m=height,
        clearance_m=curved_gap,
    )
    assert force.min_film_thickness_m == pytest.approx(expected_curved, rel=1e-9)
    assert force.min_film_thickness_m == pytest.approx(10.0e-6, abs=1e-9)


def test_reference_flat_film_reduces_to_ten_microns(geometry):
    """Zero flat approach/tilt leaves the flat film at its 10 um clearance."""

    bush = SwingBush()
    height = geometry.cylinder_height_m - 2.0 * AXIAL_CLEARANCE_M
    flat_gap, _ = film_thicknesses_m(geometry, bush)
    force = flat_slider_film(
        0.0,
        0.0,
        0.0,
        0.0,
        0.5,  # sliding speed; parallel non-approaching pad carries no load
        length_m=11.0e-3,
        height_m=height,
        clearance_m=flat_gap,
    )
    assert force.min_film_thickness_m == pytest.approx(flat_gap, rel=1e-12)
    assert force.min_film_thickness_m == pytest.approx(10.0e-6, abs=1e-9)


# ---------------------------------------------------------------------------
# Inertia bookkeeping
# ---------------------------------------------------------------------------


def test_piece_inertia_uses_parallel_axis_to_the_piece_centre(geometry):
    """Inertia about the piece centre = rescaled CM inertia + m_p d^2 (parallel axis)."""

    bush = SwingBush()
    inertia_centre = _piece_inertia_about_centre(geometry, bush, BUSH_PIECE_MASS_KG)
    props = bush_piece_mass_properties(geometry, 1.0, bush=bush)
    inertia_cg = props.inertia_about_cg_kg_m2 * (BUSH_PIECE_MASS_KG / props.mass_kg)
    # The parallel-axis term is strictly positive, so the centre inertia is larger.
    assert inertia_centre > inertia_cg
    # And matches the CM inertia plus m_p * (cg_offset - shift)^2.
    d = math.hypot(*props.cg_offset_m) - bush.piece_shift_m
    assert inertia_centre == pytest.approx(inertia_cg + BUSH_PIECE_MASS_KG * d * d, rel=1e-9)
    # Sanity magnitude (~2.7e-7 kg m^2).
    assert 1.0e-7 < inertia_centre < 1.0e-6


# ---------------------------------------------------------------------------
# Validation of inputs
# ---------------------------------------------------------------------------


def test_invalid_inputs_raise(geometry):
    with pytest.raises(ValueError):
        integrate_rotor_bush_orbit(geometry, rotor_mass_kg=0.0)
    with pytest.raises(ValueError):
        integrate_rotor_bush_orbit(geometry, piece_mass_kg=-1.0)
    with pytest.raises(ValueError):
        integrate_rotor_bush_orbit(geometry, revolutions=0)


# ---------------------------------------------------------------------------
# Steady-orbit behaviour (one slow module-scoped integration)
# ---------------------------------------------------------------------------


def test_orbit_completes_with_positive_films(orbit):
    """Every film stays open over the cycle (no metal contact)."""

    assert min(orbit.in_curved_film_m) > 0.0
    assert min(orbit.out_curved_film_m) > 0.0
    assert min(orbit.in_flat_film_m) > 0.0
    assert min(orbit.out_flat_film_m) > 0.0
    assert orbit.minimum_curved_film_m > 0.0
    assert orbit.minimum_flat_film_m > 0.0


def test_rotor_bore_stays_within_the_journal_clearance(orbit):
    """The rotor bore eccentricity ratio stays below unity (bore inside the clearance)."""

    assert 0.0 < orbit.peak_rotor_eccentricity_ratio < 1.0


def test_deviations_are_small(orbit):
    """All nine DOF stay near the prescribed reference (the perturbation premise)."""

    bush = SwingBush()
    # Rotor attitude deviation is a few mrad at most.
    assert max(abs(v) for v in orbit.rotor_attitude_deviation_rad) < 5.0e-3
    # Piece centres stay within ~15 um of their reference +/- piece_shift.
    for x in orbit.in_piece_x_m:
        assert abs(x - bush.piece_shift_m) < 20.0e-6
    for x in orbit.out_piece_x_m:
        assert abs(x + bush.piece_shift_m) < 20.0e-6


def test_orbit_is_periodic(orbit):
    """The reported revolution closes on itself (steady orbit reached)."""

    assert orbit.eccentricity_x_m[0] == pytest.approx(orbit.eccentricity_x_m[-1], abs=1.0e-6)
    assert orbit.eccentricity_y_m[0] == pytest.approx(orbit.eccentricity_y_m[-1], abs=1.0e-6)
    assert orbit.in_piece_y_m[0] == pytest.approx(orbit.in_piece_y_m[-1], abs=5.0e-5)


def test_bush_friction_is_positive_and_order_of_quasi_static(orbit, geometry, trace):
    """Dynamic bush friction is positive dissipation and O(the quasi-static estimate)."""

    quasi = friction_power_cycle_w(geometry, trace=trace)
    assert orbit.bush_friction_power_w > 0.0
    assert orbit.quasi_static_bush_friction_power_w == pytest.approx(quasi, rel=1e-9)
    # Same order of magnitude (thinner dynamic films raise it, but not wildly).
    assert 0.2 * quasi < orbit.bush_friction_power_w < 10.0 * quasi


def test_confirmed_masses_and_inertia_are_reported(orbit):
    assert orbit.piece_mass_kg == pytest.approx(BUSH_PIECE_MASS_KG)
    assert orbit.rotor_inertia_kg_m2 == pytest.approx(ROTOR_POLAR_INERTIA_KG_M2)
    assert np.isfinite(orbit.piece_inertia_kg_m2) and orbit.piece_inertia_kg_m2 > 0.0


# ---------------------------------------------------------------------------
# Rotor-cylinder sealing contact (the default full methodology, Section 4.15)
# ---------------------------------------------------------------------------


def test_seal_contact_is_engaged_and_physical(orbit):
    """The centrifugal load keeps the seal engaged; N_c is O(100s N), penetration sub-um.

    The ~6 um free-orbit bore penetration collapses to the physical Hertz deflection
    (< 3 um), and the boundary friction is real dissipation.
    """

    assert orbit.mean_seal_normal_force_n > 0.0
    assert 200.0 < orbit.mean_seal_normal_force_n < 450.0  # coupled level (rev-2 coarse)
    assert orbit.peak_seal_normal_force_n > orbit.mean_seal_normal_force_n
    assert 0.0 < orbit.max_seal_penetration_m < 3.0e-6
    assert orbit.seal_contact_friction_power_w > 0.0


def test_seal_radial_balance_and_attribution(orbit):
    """The radial force balance closes, and the load split is the attributed mechanism.

    The bush reaction is *inward* (it does not directly press the rotor out); the
    crank-pin journal, loaded by the fixed-vane constraint the bush enforces, reacts
    *outward* and is what drives the rotor onto the bore. The cylinder reacts inward.
    """

    gas, journal, bush, seal, cent = orbit.mean_radial_force_breakdown_n
    assert abs(gas + journal + bush + seal + cent) < 5.0  # steady-state, mean accel ~ 0
    assert journal > 0.0  # journal film reacts outward (the muscle)
    assert bush < 0.0  # bush reaction is inward (the messenger, not the cause)
    assert seal < 0.0  # cylinder reacts the net outward push inward
    assert 40.0 < cent < 48.0  # centrifugal m_r omega^2 e ~ 44 N, always outward
    # The seal reaction magnitude equals the reported N_c mean.
    assert -seal == pytest.approx(orbit.mean_seal_normal_force_n, rel=1e-6)


def test_gas_film_boundary_squeezes_the_flat_film(geometry, trace):
    """The 4 MPa bore gas presses the piece onto the vane (Section 4.11 gas BC).

    Referenced to the bore, the flat (vane-sealing) film runs a bore->chamber drop that
    pushes the piece toward the vane, thinning the flat film toward contact and raising
    the (Couette) bush friction sharply -- flagging the vane-bush interface for an
    EHL/contact treatment, as the rotor-cylinder seal already got.
    """

    kw = {"revolutions": 2, "samples": 40, "grid_samples": 90, "n_beta": 61, "n_s": 41,
          "trace": trace}
    hydro = integrate_rotor_bush_orbit(geometry, gas_film_boundary=False, **kw)
    gas = integrate_rotor_bush_orbit(geometry, gas_film_boundary=True, **kw)
    assert gas.minimum_flat_film_m < hydro.minimum_flat_film_m  # squeezed thinner
    assert gas.bush_friction_power_w > 3.0 * hydro.bush_friction_power_w  # much higher loss


def test_seal_contact_can_be_disabled(geometry, trace):
    """seal_contact=False isolates the bush films: no cylinder reaction, empty seal data."""

    decoupled = integrate_rotor_bush_orbit(
        geometry,
        revolutions=1,
        samples=16,
        grid_samples=90,
        n_beta=61,
        n_s=41,
        trace=trace,
        seal_contact=False,
    )
    assert decoupled.seal_normal_force_n == ()
    assert decoupled.mean_radial_force_breakdown_n == ()
    assert decoupled.mean_seal_normal_force_n == 0.0
    assert decoupled.seal_contact_friction_power_w == 0.0
    assert decoupled.bush_friction_power_w > 0.0  # the bush model still runs
