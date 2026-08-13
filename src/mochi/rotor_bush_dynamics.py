"""Nine-DOF rotor + two-piece swing-bush multibody dynamics (D2, Stage 5).

This extends the rotor lateral EOM (:mod:`mochi.rotor_dynamics`, Section 4.13) by
making the **two swing-bush pieces independent bodies** and coupling them to the rotor
through the curved (arc) and flat (slider) oil films rebuilt in the axial-uniform 1-D
form (:mod:`mochi.arc_film`, :mod:`mochi.slider_film`, Section 4.11). The shaft axis is
the cylinder axis throughout; the crank pin orbits, the shaft centreline does not.

**Perturbation about the prescribed kinematics.** The prescribed kinematics
(:func:`mochi.kinematics.prescribed_state`) is the *rigid-film limit*: the bush glued
to the vane axis and to the groove, giving the exact +/-10.37 deg rotor swing. With
compliant ~10 um films the true motion is that reference plus small (um / mrad)
deviations, so all nine DOF are modelled as **deviations about the reference** (like
the rotor centre ``e_j = O_b - O_j`` already is in Section 4.13). Zeroing every
deviation recovers the prescribed pose and the Section 4.11 10 um films exactly.

**The nine DOF** (Section 3 of the spec): rotor bore centre ``e_j = O_b - O_j``
(2), rotor attitude deviation ``dphi_r`` with ``phi_orient = phi_orient,ref +
dphi_r`` (1), and each bush piece's global centre ``O_p,k`` (2) and attitude ``phi_k``
(1) for ``k in {IN, OUT}`` -- 9 positions + 9 rates = an 18-state stiff ODE integrated
with implicit BDF, exactly as Section 4.13.

**Frames (two, deliberately).** The groove rides the rotor, so its position and the
curved-film relative geometry use the *swinging* rotor frame ``n = (cos phi_orient,
sin phi_orient)`` (the O_b->O_g direction) and ``tau = (-sin, cos)``. The vane is
*fixed* in the world (along +y), so the flat film against it is treated in the
**world** frame: the piece slides along world ``y`` and approaches the vane along
world ``x``. (Expressing the flat film in the swinging frame -- as an early draft of
the spec did -- would misstate the slide speed by the ~10 deg swing: at theta = 90 deg
the world-y slide is 0.85 m/s while the tau-component is only 0.15 m/s.)

**Curved-film clearance.** The arc film's ``clearance_m`` is the *concentric* radial
gap ``r_cut - r_b = 30 um`` (the film thickness where the piece is centred in the
groove); the reference ``piece_shift`` of 20 um toward the piece's own side then makes
the minimum curved film ``30 - 20 = 10 um`` -- matching
:func:`mochi.bush_film.film_thicknesses_m`. (This corrects the spec Section 5 slip that
listed 10 um as the arc clearance, which would drive the film straight into contact.)

See PHYSICS.md Section 4.14.
"""

from dataclasses import dataclass
from math import copysign, cos, hypot, pi, sin

import numpy as np
from scipy.integrate import solve_ivp

from mochi.arc_film import arc_film_force
from mochi.asperity_contact import AsperityParams, greenwood_tripp_pressure
from mochi.bush_film import (
    AXIAL_CLEARANCE_M,
    LUBRICANT_VISCOSITY_PA_S,
    film_thicknesses_m,
    flat_contact_length_m,
    friction_power_cycle_w,
)
from mochi.chamber_volume import SwingBush
from mochi.chambers import build_cycle_trace, seal_over_half_angle_rad
from mochi.journal_bearing import JOURNAL_CLEARANCE_M, journal_relative_speed_rad_s
from mochi.mass_properties import bush_piece_mass_properties
from mochi.ocvirk_bearing import (
    crank_pin_entrainment_speed_rad_s,
    eccentric_friction_torque_nm,
)
from mochi.rotor_cylinder import (
    ROTOR_CYLINDER_FRICTION_COEFF,
    contact_sliding_speed_m_s,
    hertz_line_contact_force_n,
)
from mochi.rotor_dynamics import _MAX_ECCENTRICITY_RATIO, ROTOR_MASS_KG, _film_force_on_rotor
from mochi.slider_film import flat_slider_film
from mochi.true_gas_force import FULL_TURN_RAD, true_gas_load

# Confirmed rotor polar inertia about its centre (user, 2026-07); paired with the
# confirmed rotor mass ``ROTOR_MASS_KG`` reused from Section 4.13.
ROTOR_POLAR_INERTIA_KG_M2 = 2.11e-4
# Confirmed single swing-bush piece mass (user, 2026-07); the raster estimate is ~6.33 g.
BUSH_PIECE_MASS_KG = 6.341e-3
# Rotor-bore / crank-side internal gas pressure (user, 2026-08): the piece is immersed in
# this, above both chambers, so it is the reference for the film gas-pressure boundary
# conditions -- only the flat (vane-sealing) film sees a bore->chamber pressure drop.
BORE_PRESSURE_PA = 4.0e6
# Mixed-lubrication asperity contact for the bush films uses the physical Greenwood-Tripp
# model (mochi.asperity_contact) -- surface parameters (sigma, beta, eta, E') rather than an
# arbitrary pressure scale -- applied to BOTH the flat (vane) pad and the curved (rotor
# groove) arc by load sharing with the hydrodynamic film. See docs/mixed_lubrication.md.
# Steel-on-steel boundary coefficient for the asperity-borne (boundary) friction.
FLAT_BOUNDARY_COEFF = 0.10

# Curved-film (rotor-groove) sealing land, from PHYSICS.md 3.3. Measured from the piece's
# own arc centre: the land runs from the recess-CHANNEL opening on one side to the rotor
# MOUTH blend on the other, giving the documented 45.0 + 41.5 = 86.5 deg land. Beyond those
# the arc faces gas, not oil, so the film -- and its gas boundary conditions -- live on the
# land: the channel end opens to the 4.0 MPa recess (= the bore reference) and the mouth end
# to the piece's chamber. The two pieces are mirror images, so the offsets swap sides.
GROOVE_LAND_MOUTH_RAD = pi * 41.5 / 180.0  # mouth (chamber) side of the land
GROOVE_LAND_CHANNEL_RAD = pi * 45.0 / 180.0  # recess-channel (bore) side of the land

# Central-difference step for the reference attitude rate/acceleration (analytic pose).
_REF_STEP_RAD = 1.0e-4
# Safety rails: keep the film solvers off metal contact during Newton iterations. The
# films are extremely stiff, so the physical state stays well inside these.
_MAX_ARC_ECC_RATIO = 0.99  # |O_p - O_g| / concentric gap
_MIN_FLAT_FILM_FRACTION = 0.02  # gap floor as a fraction of the flat clearance


@dataclass(frozen=True, slots=True)
class RotorBushOrbit:
    """Rotor + two-piece bush steady orbit over the final revolution (SI units)."""

    crank_angle_rad: tuple[float, ...]
    # Rotor bore eccentricity e_j = O_b - O_j and attitude deviation dphi_r.
    eccentricity_x_m: tuple[float, ...]
    eccentricity_y_m: tuple[float, ...]
    rotor_attitude_deviation_rad: tuple[float, ...]
    # Each piece's global centre and attitude deviation.
    in_piece_x_m: tuple[float, ...]
    in_piece_y_m: tuple[float, ...]
    in_piece_attitude_rad: tuple[float, ...]
    out_piece_x_m: tuple[float, ...]
    out_piece_y_m: tuple[float, ...]
    out_piece_attitude_rad: tuple[float, ...]
    # Minimum films per piece over the revolution (safety indicators).
    in_curved_film_m: tuple[float, ...]
    in_flat_film_m: tuple[float, ...]
    out_curved_film_m: tuple[float, ...]
    out_flat_film_m: tuple[float, ...]
    # Scalars.
    peak_rotor_eccentricity_ratio: float
    peak_curved_eccentricity_ratio: float
    minimum_curved_film_m: float
    minimum_flat_film_m: float
    bush_friction_power_w: float
    quasi_static_bush_friction_power_w: float
    curved_friction_power_w: float  # bush-rotor (curved film) component of bush_friction
    flat_friction_power_w: float  # bush-vane (flat film) component of bush_friction
    rotor_mass_kg: float
    rotor_inertia_kg_m2: float
    piece_mass_kg: float
    piece_inertia_kg_m2: float
    journal_clearance_m: float
    revolutions: int
    frequency_hz: float
    # Optional rotor-cylinder sealing contact (populated only when seal_contact=True;
    # otherwise these stay at their empty/zero defaults so the decoupled run is unchanged).
    seal_normal_force_n: tuple[float, ...] = ()
    seal_penetration_m: tuple[float, ...] = ()
    seal_contact_friction_power_w: float = 0.0
    mean_seal_normal_force_n: float = 0.0
    peak_seal_normal_force_n: float = 0.0
    max_seal_penetration_m: float = 0.0
    # Cycle-mean radial (outward, along O_b/|O_b|) force breakdown on the rotor:
    # (gas, journal film, bush reaction, seal contact, centrifugal), in newtons.
    mean_radial_force_breakdown_n: tuple[float, ...] = ()
    # Per-crank-angle raw force/moment channels for .dat export (aligned with
    # crank_angle_rad): name -> tuple of per-sample values (forces N, moments N*m). Rotor
    # forces (world x, y) per source, rotor moments about O_b, and each bush piece's film
    # force/moment. None on constructions that do not populate it.
    sample_channels: dict | None = None
    # Flat-film friction split: boundary (asperity-borne, LINEAR in FLAT_BOUNDARY_COEFF) vs
    # viscous (shear over the film), plus the cycle-mean asperity load carrying the pad. Zero
    # unless flat_mixed_lubrication is on (then boundary + viscous == flat_friction_power_w).
    flat_boundary_friction_power_w: float = 0.0
    flat_viscous_friction_power_w: float = 0.0
    mean_flat_asperity_load_n: float = 0.0


# ---------------------------------------------------------------------------
# Reference kinematics (analytic pose, cheap; differenced for rates)
# ---------------------------------------------------------------------------


def _reference_attitude(geometry, theta: float) -> tuple[float, float, float]:
    """Reference rotor orientation and its time rate/acceleration at crank angle theta.

    ``phi_orient,ref`` is :attr:`PrescribedState.rotor_orientation_rad`; the rate and
    acceleration follow from the constant drive speed, ``d/dt = omega d/dtheta`` and
    ``d2/dt2 = omega^2 d2/dtheta2``, by central differences on the analytic pose.
    """

    from mochi.kinematics import prescribed_state

    omega = geometry.angular_speed_rad_s
    h = _REF_STEP_RAD
    minus = prescribed_state(geometry, theta - h).rotor_orientation_rad
    mid = prescribed_state(geometry, theta).rotor_orientation_rad
    plus = prescribed_state(geometry, theta + h).rotor_orientation_rad
    rate = omega * (plus - minus) / (2.0 * h)
    accel = omega * omega * (plus - 2.0 * mid + minus) / (h * h)
    return mid, rate, accel


def _prescribe_gas(geometry, samples: int, trace):
    """Sample the gas force and rotor moment on a periodic grid (Section 4.13 style).

    The gas load is undefined in the small seal-over window (chambers merge); the grid
    runs over the separated range and the seam is bridged by a wrap point at 0 / 2*pi.
    Also carries the crank-pin journal entrainment speed, which is smooth everywhere.
    """

    half_angle = seal_over_half_angle_rad(geometry)
    inner = np.linspace(half_angle, FULL_TURN_RAD - half_angle, samples)
    gas_x = np.empty(samples)
    gas_y = np.empty(samples)
    moment = np.empty(samples)
    suction = np.empty(samples)
    compression = np.empty(samples)
    for index, angle in enumerate(inner):
        load = true_gas_load(geometry, float(angle), trace=trace)
        gas_x[index] = load.rotor_force_n[0]
        gas_y[index] = load.rotor_force_n[1]
        moment[index] = load.rotor_moment_about_center_nm
        suction[index] = load.suction_pressure_pa
        compression[index] = load.compression_pressure_pa
    bridge_x = 0.5 * (gas_x[0] + gas_x[-1])
    bridge_y = 0.5 * (gas_y[0] + gas_y[-1])
    bridge_m = 0.5 * (moment[0] + moment[-1])
    bridge_s = 0.5 * (suction[0] + suction[-1])
    bridge_c = 0.5 * (compression[0] + compression[-1])
    grid = np.concatenate(([0.0], inner, [FULL_TURN_RAD]))
    gx = np.concatenate(([bridge_x], gas_x, [bridge_x]))
    gy = np.concatenate(([bridge_y], gas_y, [bridge_y]))
    mz = np.concatenate(([bridge_m], moment, [bridge_m]))
    suc = np.concatenate(([bridge_s], suction, [bridge_s]))
    comp = np.concatenate(([bridge_c], compression, [bridge_c]))
    speed = np.array([crank_pin_entrainment_speed_rad_s(geometry, float(a)) for a in grid])
    return grid, gx, gy, mz, speed, suc, comp


def _piece_inertia_about_centre(geometry, bush, piece_mass_kg: float) -> float:
    """Polar inertia of one bush piece about the piece centre O_p (the arc centre).

    Uses the raster CM inertia (:func:`mochi.mass_properties.bush_piece_mass_properties`)
    rescaled to the confirmed piece mass, then the parallel-axis shift from the CM to
    the piece centre. The two pieces are mirror images, so one side suffices.
    """

    props = bush_piece_mass_properties(geometry, 1.0, bush=bush)
    mass_scale = piece_mass_kg / props.mass_kg
    inertia_cg = props.inertia_about_cg_kg_m2 * mass_scale
    # CM offset from the groove centre; the piece centre sits piece_shift toward the
    # same (own) side, so the CM->centre distance is their difference.
    cg_from_groove = hypot(*props.cg_offset_m)
    cg_to_centre = abs(cg_from_groove - bush.piece_shift_m)
    return inertia_cg + piece_mass_kg * cg_to_centre * cg_to_centre


# ---------------------------------------------------------------------------
# Coupling helpers (states -> film inputs) with contact safety rails
# ---------------------------------------------------------------------------


def _curved_film_on_piece(
    ecc_x,
    ecc_y,
    ecc_dot_x,
    ecc_dot_y,
    entrainment,
    arc_center,
    half_arc,
    radius,
    height,
    gap,
    n_beta,
    shear_speed=0.0,
    pressure_start_pa=0.0,
    pressure_end_pa=0.0,
    cavitation_pressure_pa=0.0,
):
    """Arc-film reaction on a piece, clamping the eccentricity off metal contact."""

    magnitude = hypot(ecc_x, ecc_y)
    limit = _MAX_ARC_ECC_RATIO * gap
    if magnitude > limit and magnitude > 0.0:
        scale = limit / magnitude
        ecc_x, ecc_y = ecc_x * scale, ecc_y * scale
    return arc_film_force(
        ecc_x,
        ecc_y,
        ecc_dot_x,
        ecc_dot_y,
        entrainment,
        arc_center_rad=arc_center,
        arc_half_span_rad=half_arc,
        radius_m=radius,
        length_m=height,
        clearance_m=gap,
        n_beta=n_beta,
        shear_speed_rad_s=shear_speed,
        pressure_start_pa=pressure_start_pa,
        pressure_end_pa=pressure_end_pa,
        cavitation_pressure_pa=cavitation_pressure_pa,
    )


def _flat_film_on_piece(
    approach,
    tilt,
    approach_rate,
    tilt_rate,
    slide_speed,
    length,
    height,
    clearance,
    n_s,
    *,
    pressure_start_pa=0.0,
    pressure_end_pa=0.0,
    cavitation_pressure_pa=0.0,
):
    """Flat-slider reaction on a piece, clamping approach/tilt off metal contact."""

    floor = _MIN_FLAT_FILM_FRACTION * clearance
    # Keep the central gap positive.
    approach = min(approach, clearance - floor)
    # Keep the tilted ends off the wall: |tilt| * (length / 2) <= central gap - floor.
    central = clearance - approach
    max_tilt = max(central - floor, 0.0) / (0.5 * length)
    tilt = max(-max_tilt, min(tilt, max_tilt))
    return flat_slider_film(
        approach,
        tilt,
        approach_rate,
        tilt_rate,
        slide_speed,
        length_m=length,
        height_m=height,
        clearance_m=clearance,
        n_s=n_s,
        pressure_start_pa=pressure_start_pa,
        pressure_end_pa=pressure_end_pa,
        cavitation_pressure_pa=cavitation_pressure_pa,
    )


def _flat_asperity_contact(approach, tilt, clearance, asperity, length, height, n=31):
    """Distributed Greenwood-Tripp asperity contact on the flat (vane) pad.

    Local gap ``h(s) = (clearance - approach) + s*tilt`` (the same convention as
    :func:`mochi.slider_film.flat_slider_film`); the physical asperity pressure
    (:func:`mochi.asperity_contact.greenwood_tripp_pressure`) is integrated over the pad for
    the normal force (off the vane) and the restoring moment about the pad centre.
    Distributing it lets the contact resist the piece *tilt*, not just the approach.
    """

    s = np.linspace(-0.5 * length, 0.5 * length, n)
    h = (clearance - approach) + s * tilt  # may be <= 0 (deep contact); GT handles it
    p = greenwood_tripp_pressure(h, asperity)
    ds = length / (n - 1)
    force = height * float(np.trapezoid(p, dx=ds))
    moment = height * float(np.trapezoid(p * s, dx=ds))
    return force, moment


def _curved_asperity_contact(
    ecc_x, ecc_y, arc_center, half_arc, radius, height, gap, asperity, n_beta
):
    """Distributed Greenwood-Tripp asperity contact on the curved (rotor-groove) arc.

    Local gap ``h(beta) = gap - (ecc_x cos beta + ecc_y sin beta)`` from the *real*
    (unclamped) eccentricity; the asperity pressure acts radially inward on the piece over
    ``dA = r H dbeta``, giving a contact reaction that resists the eccentricity by load
    sharing with the arc-film hydrodynamic force (replacing the hard eccentricity clamp).
    Radial pressure passes through the piece centre O_p, so it adds no moment about O_p.
    """

    beta = np.linspace(arc_center - half_arc, arc_center + half_arc, n_beta)
    cos_b, sin_b = np.cos(beta), np.sin(beta)
    h = gap - (ecc_x * cos_b + ecc_y * sin_b)  # real gap, may be <= 0
    p = greenwood_tripp_pressure(h, asperity)
    d_beta = 2.0 * half_arc / (n_beta - 1)
    force_x = -radius * height * float(np.trapezoid(p * cos_b, dx=d_beta))
    force_y = -radius * height * float(np.trapezoid(p * sin_b, dx=d_beta))
    return force_x, force_y


# ---------------------------------------------------------------------------
# The 18-state right-hand side and the integrator
# ---------------------------------------------------------------------------


def integrate_rotor_bush_orbit(
    geometry,
    *,
    rotor_mass_kg: float = ROTOR_MASS_KG,
    rotor_inertia_kg_m2: float = ROTOR_POLAR_INERTIA_KG_M2,
    piece_mass_kg: float = BUSH_PIECE_MASS_KG,
    piece_inertia_kg_m2: float | None = None,
    journal_clearance_m: float = JOURNAL_CLEARANCE_M,
    bush: SwingBush | None = None,
    revolutions: int = 4,
    samples: int = 360,
    grid_samples: int = 360,
    n_beta: int = 181,
    n_s: int = 101,
    trace=None,
    rtol: float = 1.0e-5,
    atol=None,
    velocity_atol: float = 1.0e-6,
    seal_contact: bool = True,
    seal_friction_coefficient: float = ROTOR_CYLINDER_FRICTION_COEFF,
    contact_smoothing_m: float = 3.0e-8,
    gas_film_boundary: bool = False,
    bore_pressure_pa: float = BORE_PRESSURE_PA,
    curved_shear_moment: bool = False,
    flat_mixed_lubrication: bool = False,
    roughness_m: float = 0.53e-6,
    curved_asperity: bool = True,
    asperity_params: AsperityParams | None = None,
) -> RotorBushOrbit:
    """Integrate the 9-DOF rotor + two-piece bush EOM to its steady periodic orbit.

    Runs ``revolutions`` crank turns from the prescribed reference (all deviations
    zero, piece velocities set to the reference groove velocity so the stiff squeeze
    films start unloaded) and returns the final revolution sampled at ``samples``
    points. Uses the implicit stiff integrator (BDF) -- the squeeze films are very
    stiff. ``n_beta`` / ``n_s`` set the arc / slider line-grid resolution.

    ``seal_contact`` (default **True**) adds the rotor-cylinder sealing-line contact
    (compliant Hertz line contact + boundary friction, Section 4.15) to the rotor
    lateral and attitude EOM, so the crank-pin journal, the swing bush, and the cylinder
    contact share the load by their true stiffnesses. This is the **full methodology**:
    the seal engagement is driven by the centrifugal ``m_r omega^2 e``, so it scales with
    ``omega^2`` -- any operating-speed sweep MUST keep it on or the (speed-dependent)
    rotor-cylinder loss is missing entirely. Because the bush ties the rotor to the fixed
    vane, the free-floating reduced model (:func:`mochi.rotor_cylinder.integrate_sealing_contact`)
    under-predicts ``N_c`` ~2x; here the vane constraint loads the journal, which presses
    the rotor ~1 um deeper onto the bore. Set ``seal_contact=False`` only to isolate the
    bush films (no cylinder reaction); the seal diagnostics then stay empty/zero.

    ``gas_film_boundary`` (opt-in, default **False**) imposes the chamber/crank gas
    pressures as the bush-film boundary conditions (Section 4.11): the piece is immersed in
    the bore gas ``bore_pressure_pa`` (above both chambers), so the flat (vane-sealing) film
    runs a bore->chamber pressure drop (outer end at the piece's chamber -- IN=suction,
    OUT=compression -- inner end at the bore), adding the Poiseuille gas bias and its
    throughflow leakage. Referenced to the bore pressure, the uniform immersion cancels by
    the divergence theorem, so the flat-film force *is* the net gas load on the piece; the
    curved film stays in the bore region (no drop). **It is a documented finding, not the
    default**: the bore gas drives the flat film to metal contact (bush friction ×8,
    contact-clamp-limited), so it awaits a vane-bush EHL/contact treatment (§4.11); the
    default runs the validated pure-hydrodynamic films.
    """

    geometry.validate()
    bush = SwingBush() if bush is None else bush
    bush.validate(geometry)
    if not (
        rotor_mass_kg > 0.0
        and rotor_inertia_kg_m2 > 0.0
        and piece_mass_kg > 0.0
        and journal_clearance_m > 0.0
    ):
        raise ValueError("Masses, inertia, and clearance must be positive.")
    if revolutions < 1 or samples < 8 or grid_samples < 16:
        raise ValueError("Need revolutions >= 1, samples >= 8, grid_samples >= 16.")

    trace = build_cycle_trace(geometry) if trace is None else trace
    omega = geometry.angular_speed_rad_s
    throw = geometry.eccentricity_m
    lever = geometry.cutout_offset_m  # groove-centre lever from the rotor bore, l_g
    radius = bush.piece_outer_radius_m
    half_arc = bush.half_arc_rad()
    height = geometry.cylinder_height_m - 2.0 * AXIAL_CLEARANCE_M
    flat_gap, _ = film_thicknesses_m(geometry, bush)  # 10 um flat clearance c_f
    curved_gap = geometry.cutout_radius_m - radius  # 30 um concentric radial gap
    shift = bush.piece_shift_m
    cyl_height = geometry.cylinder_height_m  # Hertz line-contact length (full rotor axial)
    rotor_radius = geometry.rotor_radius_m  # lever for the seal-friction moment about O_b
    m_r, i_r, m_p = rotor_mass_kg, rotor_inertia_kg_m2, piece_mass_kg
    i_p = (
        _piece_inertia_about_centre(geometry, bush, piece_mass_kg)
        if piece_inertia_kg_m2 is None
        else piece_inertia_kg_m2
    )
    # Greenwood-Tripp asperity parameters for the mixed-lubrication films (both bush films).
    asperity = (
        asperity_params if asperity_params is not None else AsperityParams(roughness_m=roughness_m)
    )
    asperity.validate()

    grid, gas_x, gas_y, gas_m, speed, suc_p, comp_p = _prescribe_gas(geometry, grid_samples, trace)

    def gas_at(theta):
        wrapped = theta % FULL_TURN_RAD
        return (
            float(np.interp(wrapped, grid, gas_x)),
            float(np.interp(wrapped, grid, gas_y)),
            float(np.interp(wrapped, grid, gas_m)),
        )

    def journal_entrainment_at(theta):
        return float(np.interp(theta % FULL_TURN_RAD, grid, speed))

    def chamber_pressures_at(theta):
        """(suction, compression) absolute gas pressures bounding the vane at ``theta``."""
        wrapped = theta % FULL_TURN_RAD
        return (
            float(np.interp(wrapped, grid, suc_p)),
            float(np.interp(wrapped, grid, comp_p)),
        )

    # side, and index offsets of (x, y, phi) in the position block
    pieces = ((1.0, 3, 4, 5), (-1.0, 6, 7, 8))

    def evaluate(theta, state, *, want_films=False):
        """Return the 9 accelerations (and optionally per-piece film diagnostics)."""

        e_jx, e_jy, dphi_r = state[0], state[1], state[2]
        ve_jx, ve_jy, dphi_r_dot = state[9], state[10], state[11]

        phi_ref, phi_dot_ref, phi_ddot_ref = _reference_attitude(geometry, theta)
        phi = phi_ref + dphi_r
        phi_dot = phi_dot_ref + dphi_r_dot
        n_hat = (cos(phi), sin(phi))
        tau_hat = (-sin(phi), cos(phi))

        pin = (throw * sin(theta), throw * cos(theta))
        pin_dot = (throw * omega * cos(theta), -throw * omega * sin(theta))
        bore = (pin[0] + e_jx, pin[1] + e_jy)
        bore_dot = (pin_dot[0] + ve_jx, pin_dot[1] + ve_jy)
        groove = (bore[0] + lever * n_hat[0], bore[1] + lever * n_hat[1])
        groove_dot = (
            bore_dot[0] + lever * phi_dot * tau_hat[0],
            bore_dot[1] + lever * phi_dot * tau_hat[1],
        )

        gx, gy, gm = gas_at(theta)
        entrain_j = journal_entrainment_at(theta)
        f_jx, f_jy = _film_force_on_rotor(e_jx, e_jy, ve_jx, ve_jy, entrain_j, journal_clearance_m)

        flat_length = flat_contact_length_m(geometry, theta, bush)

        # Chamber gas pressures at the vane (for the flat-film gas boundary condition).
        # The piece is immersed in the bore gas (bore_pressure_pa) so the film pressures
        # are referenced to it; only the flat (vane-sealing) film sees a bore->chamber
        # drop, with the outer end at the piece's own chamber (IN=+x=suction,
        # OUT=-x=compression) and the inner end at the bore. Gauge-relative-to-bore, so
        # the uniform bore pressure drops out of the net piece force (divergence theorem).
        suction_pa, compression_pa = chamber_pressures_at(theta)

        curved_reaction_x = 0.0
        curved_reaction_y = 0.0
        rotor_curved_moment = 0.0
        piece_accels = []
        films = []
        piece_detail = []  # per-piece (curved F, flat N, moments) for the sample diagnostics
        for side, ix, iy, ip in pieces:
            x_k, y_k, phi_k = state[ix], state[iy], state[ip]
            vx_k, vy_k, vphi_k = state[ix + 9], state[iy + 9], state[ip + 9]

            # Curved film (piece OD in the rotor groove); ecc = O_p - O_g (world).
            ecc_x, ecc_y = x_k - groove[0], y_k - groove[1]
            ecc_dot_x, ecc_dot_y = vx_k - groove_dot[0], vy_k - groove_dot[1]
            arc_center = (0.0 if side > 0.0 else pi) + phi_k
            entrainment_c = 0.5 * (phi_dot + vphi_k)  # mean groove+piece spin about O_g
            # Groove-relative-to-piece rotation shears the curved film (rotor swing drag).
            shear_c = (phi_dot - vphi_k) if curved_shear_moment else 0.0

            # Curved-film gas boundary (opt-in with gas_film_boundary): the film lives on the
            # sealing LAND, which runs from the recess-channel opening (4.0 MPa = the bore
            # reference, gauge 0) to the mouth blend (the piece's own chamber). The land is
            # asymmetric about the arc centre and mirrors between the two pieces, so it is
            # expressed as a shifted centre + half-span. Beyond the land the arc faces gas,
            # not oil, so restricting the film to the land is itself the physical geometry.
            arc_c, arc_half = arc_center, half_arc
            gas_start_c = gas_end_c = cav_c = 0.0
            if gas_film_boundary:
                chamber_c = suction_pa if side > 0.0 else compression_pa
                chamber_gauge_c = chamber_c - bore_pressure_pa
                cav_c = min(0.0, chamber_gauge_c)  # same floor rule as the flat film
                arc_half = 0.5 * (GROOVE_LAND_MOUTH_RAD + GROOVE_LAND_CHANNEL_RAD)
                offset = 0.5 * (GROOVE_LAND_MOUTH_RAD - GROOVE_LAND_CHANNEL_RAD)
                if side > 0.0:  # IN: channel at -45 deg, mouth at +41.5 deg
                    arc_c = arc_center + offset
                    gas_start_c, gas_end_c = 0.0, chamber_gauge_c
                else:  # OUT: mirrored -- mouth at -41.5 deg, channel at +45 deg
                    arc_c = arc_center - offset
                    gas_start_c, gas_end_c = chamber_gauge_c, 0.0

            arc = _curved_film_on_piece(
                ecc_x,
                ecc_y,
                ecc_dot_x,
                ecc_dot_y,
                entrainment_c,
                arc_c,
                arc_half,
                radius,
                height,
                curved_gap,
                n_beta,
                shear_speed=shear_c,
                pressure_start_pa=gas_start_c,
                pressure_end_pa=gas_end_c,
                cavitation_pressure_pa=cav_c,
            )
            fc_x, fc_y = arc.force_x_n, arc.force_y_n  # on the piece

            # Flat film (piece flat on the fixed vane) in the world frame: approach
            # along world x (toward the vane), slide along world y. The pad's ``s`` runs
            # from -L/2 (inner, radially toward O) to +L/2 (outer, toward the chamber), so
            # pressure_start = bore (gauge 0), pressure_end = chamber - bore (< 0).
            approach = shift - side * x_k  # 0 at the reference x_k = side * shift
            approach_rate = -side * vx_k
            outer_gauge = 0.0
            cav_gauge = 0.0
            if gas_film_boundary:
                chamber_pa = suction_pa if side > 0.0 else compression_pa
                outer_gauge = chamber_pa - bore_pressure_pa
                # Cavitation floor = the LOWEST pressure the film is connected to (here the
                # chamber). Refrigerant-saturated oil evolves gas as soon as it drops below
                # the local saturation pressure, so a film cannot sustain a pressure below
                # its own feed; this is the same rule the journal film already follows (both
                # its ends sit at the 4.0 MPa ambient, so its Gumbel floor IS that ambient --
                # PHYSICS.md 4.10). Using absolute zero instead would let the diverging
                # half-stroke of this reciprocating pad hold a multi-MPa suction lobe.
                cav_gauge = min(0.0, outer_gauge)
            slider = _flat_film_on_piece(
                approach,
                phi_k,
                approach_rate,
                vphi_k,
                vy_k,
                flat_length,
                height,
                flat_gap,
                n_s,
                pressure_start_pa=0.0,
                pressure_end_pa=outer_gauge,
                cavitation_pressure_pa=cav_gauge,
            )
            f_flat = side * slider.normal_force_n  # off the vane, along +side * x_hat

            # Mixed lubrication (both films): Greenwood-Tripp asperity contact carries part of
            # the load by SHARING with the hydrodynamic film, from the *real* (unclamped)
            # separation -- a compliant contact replacing the hard clamps. Flat pad: acts off
            # the vane (+ tilt-restoring moment). Curved arc: acts radially inward on the piece
            # (resists the eccentricity), with an equal-opposite reaction on the rotor groove.
            asp_force = asp_fx = asp_moment = 0.0
            asp_c_x = asp_c_y = 0.0
            if flat_mixed_lubrication:
                asp_force, asp_moment = _flat_asperity_contact(
                    approach, phi_k, flat_gap, asperity, flat_length, height
                )
                asp_fx = side * asp_force
                if curved_asperity:
                    asp_c_x, asp_c_y = _curved_asperity_contact(
                        ecc_x, ecc_y, arc_center, half_arc, radius, height, curved_gap,
                        asperity, n_beta,
                    )

            # Piece equations of motion (absolute world frame). The gas-pressure boundary
            # condition on the flat film carries the net gas load on the piece (referenced
            # to the bore pressure, so the immersing bore gas cancels by the divergence
            # theorem); the curved film stays in the bore region (no gas drop).
            piece_accels.append(
                (
                    (fc_x + asp_c_x + f_flat + asp_fx) / m_p,
                    (fc_y + asp_c_y) / m_p,
                    (slider.moment_nm + arc.shear_moment_nm + asp_moment) / i_p,
                )
            )

            # Rotor reaction: -(F_c + F_asp,curved) on the groove at O_g, lever (O_g - O_b) =
            # l_g * n. The curved shear torque on the piece has an equal-opposite rotor reaction.
            curved_reaction_x -= fc_x + asp_c_x
            curved_reaction_y -= fc_y + asp_c_y
            react_x, react_y = -(fc_x + asp_c_x), -(fc_y + asp_c_y)
            rotor_curved_moment += (lever * n_hat[0]) * react_y - (lever * n_hat[1]) * react_x
            rotor_curved_moment -= arc.shear_moment_nm

            if want_films:
                films.append((arc.min_film_thickness_m, slider.min_film_thickness_m))
                piece_detail.append(
                    (
                        fc_x,  # curved-film force on the piece, x (world)
                        fc_y,  # curved-film force on the piece, y (world)
                        slider.normal_force_n,  # flat-film normal force (off the vane)
                        slider.moment_nm,  # flat-film (slider) moment about O_p
                        arc.shear_moment_nm,  # curved-film Couette shear moment about O_p
                        asp_force,  # asperity contact normal force (mixed lubrication)
                        asp_moment,  # asperity contact moment about O_p
                    )
                )

        # Optional rotor-cylinder sealing contact (compliant Hertz line contact + boundary
        # friction), ported from mochi.rotor_cylinder. The OD touches the bore where the
        # rotor is farthest out: penetration = |O_b| - throw (R_cyl - R_r = throw), normal
        # along u = O_b/|O_b| (radially inward reaction), friction tangential at O_b + R_r u.
        seal_fx = seal_fy = seal_moment = seal_load = 0.0
        if seal_contact:
            bore_r = hypot(bore[0], bore[1])
            u_sx, u_sy = bore[0] / bore_r, bore[1] / bore_r
            overlap = bore_r - throw
            softened = 0.5 * (overlap + hypot(overlap, contact_smoothing_m))  # smooth max(.,0)
            seal_load = hertz_line_contact_force_n(softened, cyl_height)
            tangent = -copysign(1.0, contact_sliding_speed_m_s(geometry, theta))
            fric = seal_friction_coefficient * seal_load * tangent
            seal_fx = -seal_load * u_sx - fric * (-u_sy)
            seal_fy = -seal_load * u_sy - fric * u_sx
            seal_moment = -rotor_radius * fric  # friction lever R_r about O_b (normal: 0)

        # Rotor lateral EOM: m_r e_j'' = F_gas + F_journal + sum(-F_c) + F_seal + m_r omega^2 O_j.
        acc_e_jx = (gx + f_jx + curved_reaction_x + seal_fx + m_r * omega * omega * pin[0]) / m_r
        acc_e_jy = (gy + f_jy + curved_reaction_y + seal_fy + m_r * omega * omega * pin[1]) / m_r

        # Rotor attitude EOM about O_b: I_r phi_orient'' = M_gas + M_fric + curved + seal,
        # with phi_orient'' = phi_ddot_ref + dphi_r'' -> subtract the reference inertia.
        eps_j = min(hypot(e_jx, e_jy) / journal_clearance_m, _MAX_ECCENTRICITY_RATIO)
        omega_rel = journal_relative_speed_rad_s(geometry, theta)
        m_fric = -copysign(eccentric_friction_torque_nm(eps_j, omega_rel), phi_dot)
        acc_dphi_r = (gm + m_fric + rotor_curved_moment + seal_moment) / i_r - phi_ddot_ref

        accels = [
            acc_e_jx,
            acc_e_jy,
            acc_dphi_r,
            piece_accels[0][0],
            piece_accels[0][1],
            piece_accels[0][2],
            piece_accels[1][0],
            piece_accels[1][1],
            piece_accels[1][2],
        ]
        if want_films:
            # Force breakdown for the seal-coupling diagnostics: gas, journal film, bush
            # reaction sum(-F_c), seal contact (vector + N_c magnitude), centrifugal.
            forces = (
                gx,
                gy,
                f_jx,
                f_jy,
                curved_reaction_x,
                curved_reaction_y,
                seal_fx,
                seal_fy,
                seal_load,
                m_r * omega * omega * pin[0],
                m_r * omega * omega * pin[1],
            )
            # Rotor moments about O_b: gas, eccentric journal friction, curved bush
            # reaction, seal-contact friction (all N*m, same sign convention as the EOM).
            moments = (gm, m_fric, rotor_curved_moment, seal_moment)
            return accels, films, forces, moments, piece_detail
        return accels

    def rhs(t, state):
        theta = omega * t
        accels = evaluate(theta, state)
        return list(state[9:]) + accels

    # --- Initial condition: the prescribed reference at the start angle ---
    from mochi.kinematics import prescribed_state

    start_theta = grid[1]  # first separated angle after the seal-over seam
    ref = prescribed_state(geometry, start_theta)
    groove_y0 = ref.cutout_center_m[1]
    phi0, ref_rate, _ = _reference_attitude(geometry, start_theta)
    # Reference groove velocity O_g' = O_j' + l_g phi_ref' tau, computed with the SAME
    # expression the RHS uses; the pieces start matching it so the stiff squeeze films
    # begin unloaded (ecc_dot = O_p' - O_g' = 0 at t_start).
    groove_vx0 = throw * omega * cos(start_theta) + lever * ref_rate * (-sin(phi0))
    groove_vy0 = -throw * omega * sin(start_theta) + lever * ref_rate * cos(phi0)

    # Warm-start the flat-film mixed contact: without it the gas load slowly drains the
    # 10 um film to its near-contact equilibrium over ~15 revolutions (a squeeze-limited
    # transient). The Greenwood-Tripp contact stiffens steeply near ~sigma, so a gas-loaded
    # pad settles at h ~ 1-2 sigma; start each loaded piece at 1.5 sigma to skip the slow
    # transient (the integration refines the exact value). Reference (unloaded) IC otherwise.
    x_in0, x_out0 = shift, -shift
    if gas_film_boundary and flat_mixed_lubrication:
        suc0, comp0 = chamber_pressures_at(start_theta)
        h_eq = min(max(1.5 * asperity.roughness_m, _MIN_FLAT_FILM_FRACTION * flat_gap), flat_gap)

        def _approach_eq(chamber_pa):
            dp = abs(chamber_pa - bore_pressure_pa)
            return 0.0 if dp <= 0.0 else (flat_gap - h_eq)

        x_in0 = shift - _approach_eq(suc0)
        x_out0 = -(shift - _approach_eq(comp0))

    state0 = [
        0.0,  # e_jx
        0.0,  # e_jy
        0.0,  # dphi_r
        x_in0,  # x_IN (warm-started under the gas load when flat mixed)
        groove_y0,  # y_IN
        0.0,  # phi_IN
        x_out0,  # x_OUT
        groove_y0,  # y_OUT
        0.0,  # phi_OUT
        0.0,  # e_jx'
        0.0,  # e_jy'
        0.0,  # dphi_r'
        groove_vx0,  # x_IN'
        groove_vy0,  # y_IN'
        0.0,  # phi_IN'
        groove_vx0,  # x_OUT'
        groove_vy0,  # y_OUT'
        0.0,  # phi_OUT'
    ]

    t_start = start_theta / omega
    period = 1.0 / geometry.frequency_hz
    t_end = t_start + revolutions * period
    t_eval = np.linspace(t_end - period, t_end, samples)
    t_eval[0] = max(t_eval[0], t_start)  # guard FP underflow at revolutions == 1

    # The squeeze films make the piece velocity modes extremely stiff (a fast,
    # overdamped ~3e-8 s eigenvalue). A uniformly tight ``atol`` pins BDF to that
    # timescale forever; a *loose* absolute tolerance on the 9 velocity states (they
    # are only a means to the positions) lets BDF step over the decayed fast mode,
    # while the positions -- the eccentricities and films we report -- stay tight.
    if atol is None:
        atol = np.full(18, 1.0e-10)
        atol[0:2] = 1.0e-11  # rotor bore eccentricity (sub-nm)
        atol[9:18] = velocity_atol  # loose on all 9 velocities
    solution = solve_ivp(
        rhs,
        (t_start, t_end),
        state0,
        method="BDF",
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
    )
    if not solution.success:
        raise RuntimeError(f"Rotor-bush orbit integration failed: {solution.message}")

    y = solution.y
    angles = (omega * t_eval) % FULL_TURN_RAD

    # Re-evaluate the films at each output sample for the min-film / friction diagnostics.
    in_curved = np.empty(samples)
    in_flat = np.empty(samples)
    out_curved = np.empty(samples)
    out_flat = np.empty(samples)
    curved_ecc_ratio = np.empty(samples)
    friction = np.empty(samples)
    curved_friction = np.empty(samples)
    flat_friction = np.empty(samples)
    # Flat-film friction split into its boundary (asperity-borne, linear in mu_b) and viscous
    # (shear over the film) parts, plus the asperity load itself -- the calibration handles.
    flat_boundary_friction = np.zeros(samples)
    flat_viscous_friction = np.zeros(samples)
    mean_asperity_load = np.zeros(samples)
    seal_normal = np.zeros(samples)
    seal_pen = np.zeros(samples)
    seal_fric = np.zeros(samples)
    # Radial (along u = O_b/|O_b|) force breakdown for the seal-coupling attribution:
    # columns = gas, journal film, bush reaction, seal contact, centrifugal.
    radial_breakdown = np.zeros((samples, 5))
    # Per-crank-angle force/moment channels for raw-data (.dat) export. Rotor forces (world
    # x, y) per source, rotor moments about O_b, and each bush piece's film force/moment.
    _channel_names = (
        "gas_force_x_n", "gas_force_y_n",
        "journal_force_x_n", "journal_force_y_n",
        "bush_reaction_x_n", "bush_reaction_y_n",
        "seal_force_x_n", "seal_force_y_n",
        "centrifugal_x_n", "centrifugal_y_n",
        "gas_moment_nm", "journal_friction_moment_nm",
        "bush_curved_moment_nm", "seal_moment_nm",
        "in_curved_force_x_n", "in_curved_force_y_n", "in_flat_normal_n",
        "in_flat_moment_nm", "in_curved_shear_moment_nm",
        "in_asperity_force_n", "in_asperity_moment_nm",
        "out_curved_force_x_n", "out_curved_force_y_n", "out_flat_normal_n",
        "out_flat_moment_nm", "out_curved_shear_moment_nm",
        "out_asperity_force_n", "out_asperity_moment_nm",
    )
    channels = {name: np.zeros(samples) for name in _channel_names}
    mu = LUBRICANT_VISCOSITY_PA_S
    arc_length = 2.0 * radius * half_arc
    for j in range(samples):
        theta = float(omega * t_eval[j])
        state = y[:, j]
        _, films, forces, moments, piece_detail = evaluate(theta, state, want_films=True)
        in_curved[j], in_flat[j] = films[0]
        out_curved[j], out_flat[j] = films[1]

        phi_ref, phi_dot_ref, _ = _reference_attitude(geometry, theta)
        phi = phi_ref + state[2]
        phi_dot = phi_dot_ref + state[11]
        bore = (throw * sin(theta) + state[0], throw * cos(theta) + state[1])
        groove = (bore[0] + lever * cos(phi), bore[1] + lever * sin(phi))
        ecc_in = hypot(state[3] - groove[0], state[4] - groove[1])
        ecc_out = hypot(state[6] - groove[0], state[7] - groove[1])
        curved_ecc_ratio[j] = max(ecc_in, ecc_out) / curved_gap

        if seal_contact:
            overlap = hypot(bore[0], bore[1]) - throw
            softened = 0.5 * (overlap + hypot(overlap, contact_smoothing_m))
            seal_normal[j] = hertz_line_contact_force_n(softened, cyl_height)
            seal_pen[j] = overlap
            seal_fric[j] = (
                seal_friction_coefficient
                * seal_normal[j]
                * abs(contact_sliding_speed_m_s(geometry, theta))
            )

        # Radial force breakdown (project each rotor force onto u = O_b/|O_b|, outward).
        bore_r = hypot(bore[0], bore[1])
        u_rx, u_ry = bore[0] / bore_r, bore[1] / bore_r
        gx_f, gy_f, fjx_f, fjy_f, brx_f, bry_f, sfx_f, sfy_f, _sl, cfx_f, cfy_f = forces
        radial_breakdown[j, 0] = gx_f * u_rx + gy_f * u_ry  # gas
        radial_breakdown[j, 1] = fjx_f * u_rx + fjy_f * u_ry  # journal film
        radial_breakdown[j, 2] = brx_f * u_rx + bry_f * u_ry  # bush reaction sum(-F_c)
        radial_breakdown[j, 3] = sfx_f * u_rx + sfy_f * u_ry  # seal contact
        radial_breakdown[j, 4] = cfx_f * u_rx + cfy_f * u_ry  # centrifugal

        # Raw per-angle force/moment channels (world frame) for the .dat export.
        gm_f, mfric_f, mcurved_f, mseal_f = moments
        (in_fcx, in_fcy, in_fn, in_mfl, in_msh, in_fasp, in_masp) = piece_detail[0]
        (out_fcx, out_fcy, out_fn, out_mfl, out_msh, out_fasp, out_masp) = piece_detail[1]
        channels["gas_force_x_n"][j] = gx_f
        channels["gas_force_y_n"][j] = gy_f
        channels["journal_force_x_n"][j] = fjx_f
        channels["journal_force_y_n"][j] = fjy_f
        channels["bush_reaction_x_n"][j] = brx_f
        channels["bush_reaction_y_n"][j] = bry_f
        channels["seal_force_x_n"][j] = sfx_f
        channels["seal_force_y_n"][j] = sfy_f
        channels["centrifugal_x_n"][j] = cfx_f
        channels["centrifugal_y_n"][j] = cfy_f
        channels["gas_moment_nm"][j] = gm_f
        channels["journal_friction_moment_nm"][j] = mfric_f
        channels["bush_curved_moment_nm"][j] = mcurved_f
        channels["seal_moment_nm"][j] = mseal_f
        channels["in_curved_force_x_n"][j] = in_fcx
        channels["in_curved_force_y_n"][j] = in_fcy
        channels["in_flat_normal_n"][j] = in_fn
        channels["in_flat_moment_nm"][j] = in_mfl
        channels["in_curved_shear_moment_nm"][j] = in_msh
        channels["in_asperity_force_n"][j] = in_fasp
        channels["in_asperity_moment_nm"][j] = in_masp
        channels["out_curved_force_x_n"][j] = out_fcx
        channels["out_curved_force_y_n"][j] = out_fcy
        channels["out_flat_normal_n"][j] = out_fn
        channels["out_flat_moment_nm"][j] = out_mfl
        channels["out_curved_shear_moment_nm"][j] = out_msh
        channels["out_asperity_force_n"][j] = out_fasp
        channels["out_asperity_moment_nm"][j] = out_masp

        # Couette bush-film friction (comparable to bush_film.friction_power_cycle_w):
        # curved sliding = (rotor swing - piece spin) * r_b; flat sliding = piece world-y.
        flat_length = flat_contact_length_m(geometry, theta, bush)
        curved_p = flat_p = 0.0
        flat_boundary_p = flat_viscous_p = 0.0
        asperity_load = 0.0
        for side_f, ix, iy, ip, cfilm, ffilm in (
            (1.0, 3, 4, 5, in_curved[j], in_flat[j]),
            (-1.0, 6, 7, 8, out_curved[j], out_flat[j]),
        ):
            u_curved = abs(phi_dot - state[ip + 9]) * radius
            curved_p += mu * u_curved * u_curved / cfilm * arc_length * height
            u_flat = abs(state[iy + 9])
            if flat_mixed_lubrication:
                # Mixed-film friction (both films): boundary (asperity load * mu_b) + viscous
                # shear over the film. The asperity load is the distributed compliant contact.
                # The two are tracked separately: the boundary term is LINEAR in mu_b and in
                # the asperity load, so a calibration can be solved for directly.
                approach_f = shift - side_f * state[ix]
                f_asp, _ = _flat_asperity_contact(
                    approach_f, state[ip], flat_gap, asperity, flat_length, height
                )
                boundary_term = FLAT_BOUNDARY_COEFF * f_asp * u_flat
                viscous_term = mu * u_flat * u_flat / ffilm * flat_length * height
                flat_boundary_p += boundary_term
                flat_viscous_p += viscous_term
                asperity_load += f_asp
                flat_p += boundary_term + viscous_term
                # Curved boundary friction (asperity-borne); tiny when the piece co-rotates
                # with the rotor (u_curved ~ 0) but retained for completeness.
                if curved_asperity:
                    arc_center_f = (0.0 if side_f > 0.0 else pi) + state[ip]
                    fca_x, fca_y = _curved_asperity_contact(
                        state[ix] - groove[0], state[iy] - groove[1],
                        arc_center_f, half_arc, radius, height, curved_gap, asperity, n_beta,
                    )
                    curved_p += FLAT_BOUNDARY_COEFF * hypot(fca_x, fca_y) * u_curved
            else:
                viscous_term = mu * u_flat * u_flat / ffilm * flat_length * height
                flat_viscous_p += viscous_term
                flat_p += viscous_term
        curved_friction[j] = curved_p
        flat_friction[j] = flat_p
        flat_boundary_friction[j] = flat_boundary_p
        flat_viscous_friction[j] = flat_viscous_p
        mean_asperity_load[j] = asperity_load
        friction[j] = curved_p + flat_p

    bush_friction = float(np.mean(friction))
    quasi_static = friction_power_cycle_w(geometry, trace=trace)

    peak_rotor_ecc = float(np.max(np.hypot(y[0], y[1])) / journal_clearance_m)

    sample_channels = {name: tuple(float(v) for v in arr) for name, arr in channels.items()}

    return RotorBushOrbit(
        crank_angle_rad=tuple(float(a) for a in angles),
        eccentricity_x_m=tuple(float(v) for v in y[0]),
        eccentricity_y_m=tuple(float(v) for v in y[1]),
        rotor_attitude_deviation_rad=tuple(float(v) for v in y[2]),
        in_piece_x_m=tuple(float(v) for v in y[3]),
        in_piece_y_m=tuple(float(v) for v in y[4]),
        in_piece_attitude_rad=tuple(float(v) for v in y[5]),
        out_piece_x_m=tuple(float(v) for v in y[6]),
        out_piece_y_m=tuple(float(v) for v in y[7]),
        out_piece_attitude_rad=tuple(float(v) for v in y[8]),
        in_curved_film_m=tuple(float(v) for v in in_curved),
        in_flat_film_m=tuple(float(v) for v in in_flat),
        out_curved_film_m=tuple(float(v) for v in out_curved),
        out_flat_film_m=tuple(float(v) for v in out_flat),
        peak_rotor_eccentricity_ratio=peak_rotor_ecc,
        peak_curved_eccentricity_ratio=float(np.max(curved_ecc_ratio)),
        minimum_curved_film_m=float(min(in_curved.min(), out_curved.min())),
        minimum_flat_film_m=float(min(in_flat.min(), out_flat.min())),
        bush_friction_power_w=bush_friction,
        quasi_static_bush_friction_power_w=quasi_static,
        curved_friction_power_w=float(np.mean(curved_friction)),
        flat_friction_power_w=float(np.mean(flat_friction)),
        rotor_mass_kg=m_r,
        rotor_inertia_kg_m2=i_r,
        piece_mass_kg=m_p,
        piece_inertia_kg_m2=i_p,
        journal_clearance_m=journal_clearance_m,
        revolutions=revolutions,
        frequency_hz=geometry.frequency_hz,
        seal_normal_force_n=tuple(float(v) for v in seal_normal) if seal_contact else (),
        seal_penetration_m=tuple(float(v) for v in seal_pen) if seal_contact else (),
        seal_contact_friction_power_w=float(np.mean(seal_fric)) if seal_contact else 0.0,
        mean_seal_normal_force_n=float(np.mean(seal_normal)) if seal_contact else 0.0,
        peak_seal_normal_force_n=float(np.max(seal_normal)) if seal_contact else 0.0,
        max_seal_penetration_m=float(np.max(seal_pen)) if seal_contact else 0.0,
        mean_radial_force_breakdown_n=(
            tuple(float(v) for v in radial_breakdown.mean(axis=0)) if seal_contact else ()
        ),
        sample_channels=sample_channels,
        flat_boundary_friction_power_w=float(np.mean(flat_boundary_friction)),
        flat_viscous_friction_power_w=float(np.mean(flat_viscous_friction)),
        mean_flat_asperity_load_n=float(np.mean(mean_asperity_load)),
    )
