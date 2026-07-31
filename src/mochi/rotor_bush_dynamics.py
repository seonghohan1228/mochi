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
from mochi.rotor_dynamics import _MAX_ECCENTRICITY_RATIO, ROTOR_MASS_KG, _film_force_on_rotor
from mochi.slider_film import flat_slider_film
from mochi.true_gas_force import FULL_TURN_RAD, true_gas_load

# Confirmed rotor polar inertia about its centre (user, 2026-07); paired with the
# confirmed rotor mass ``ROTOR_MASS_KG`` reused from Section 4.13.
ROTOR_POLAR_INERTIA_KG_M2 = 2.11e-4
# Confirmed single swing-bush piece mass (user, 2026-07); the raster estimate is ~6.33 g.
BUSH_PIECE_MASS_KG = 6.341e-3

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
    rotor_mass_kg: float
    rotor_inertia_kg_m2: float
    piece_mass_kg: float
    piece_inertia_kg_m2: float
    journal_clearance_m: float
    revolutions: int
    frequency_hz: float


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
    for index, angle in enumerate(inner):
        load = true_gas_load(geometry, float(angle), trace=trace)
        gas_x[index] = load.rotor_force_n[0]
        gas_y[index] = load.rotor_force_n[1]
        moment[index] = load.rotor_moment_about_center_nm
    bridge_x = 0.5 * (gas_x[0] + gas_x[-1])
    bridge_y = 0.5 * (gas_y[0] + gas_y[-1])
    bridge_m = 0.5 * (moment[0] + moment[-1])
    grid = np.concatenate(([0.0], inner, [FULL_TURN_RAD]))
    gx = np.concatenate(([bridge_x], gas_x, [bridge_x]))
    gy = np.concatenate(([bridge_y], gas_y, [bridge_y]))
    mz = np.concatenate(([bridge_m], moment, [bridge_m]))
    speed = np.array([crank_pin_entrainment_speed_rad_s(geometry, float(a)) for a in grid])
    return grid, gx, gy, mz, speed


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
    )


def _flat_film_on_piece(
    approach, tilt, approach_rate, tilt_rate, slide_speed, length, height, clearance, n_s
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
    )


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
) -> RotorBushOrbit:
    """Integrate the 9-DOF rotor + two-piece bush EOM to its steady periodic orbit.

    Runs ``revolutions`` crank turns from the prescribed reference (all deviations
    zero, piece velocities set to the reference groove velocity so the stiff squeeze
    films start unloaded) and returns the final revolution sampled at ``samples``
    points. Uses the implicit stiff integrator (BDF) -- the squeeze films are very
    stiff. ``n_beta`` / ``n_s`` set the arc / slider line-grid resolution.
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
    m_r, i_r, m_p = rotor_mass_kg, rotor_inertia_kg_m2, piece_mass_kg
    i_p = (
        _piece_inertia_about_centre(geometry, bush, piece_mass_kg)
        if piece_inertia_kg_m2 is None
        else piece_inertia_kg_m2
    )

    grid, gas_x, gas_y, gas_m, speed = _prescribe_gas(geometry, grid_samples, trace)

    def gas_at(theta):
        wrapped = theta % FULL_TURN_RAD
        return (
            float(np.interp(wrapped, grid, gas_x)),
            float(np.interp(wrapped, grid, gas_y)),
            float(np.interp(wrapped, grid, gas_m)),
        )

    def journal_entrainment_at(theta):
        return float(np.interp(theta % FULL_TURN_RAD, grid, speed))

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

        curved_reaction_x = 0.0
        curved_reaction_y = 0.0
        rotor_curved_moment = 0.0
        piece_accels = []
        films = []
        for side, ix, iy, ip in pieces:
            x_k, y_k, phi_k = state[ix], state[iy], state[ip]
            vx_k, vy_k, vphi_k = state[ix + 9], state[iy + 9], state[ip + 9]

            # Curved film (piece OD in the rotor groove); ecc = O_p - O_g (world).
            ecc_x, ecc_y = x_k - groove[0], y_k - groove[1]
            ecc_dot_x, ecc_dot_y = vx_k - groove_dot[0], vy_k - groove_dot[1]
            arc_center = (0.0 if side > 0.0 else pi) + phi_k
            entrainment_c = 0.5 * (phi_dot + vphi_k)  # mean groove+piece spin about O_g
            arc = _curved_film_on_piece(
                ecc_x,
                ecc_y,
                ecc_dot_x,
                ecc_dot_y,
                entrainment_c,
                arc_center,
                half_arc,
                radius,
                height,
                curved_gap,
                n_beta,
            )
            fc_x, fc_y = arc.force_x_n, arc.force_y_n  # on the piece

            # Flat film (piece flat on the fixed vane) in the world frame: approach
            # along world x (toward the vane), slide along world y.
            approach = shift - side * x_k  # 0 at the reference x_k = side * shift
            approach_rate = -side * vx_k
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
            )
            f_flat = side * slider.normal_force_n  # off the vane, along +side * x_hat

            # Piece equations of motion (absolute world frame; gas bias omitted).
            piece_accels.append(
                (
                    (fc_x + f_flat) / m_p,
                    fc_y / m_p,
                    slider.moment_nm / i_p,
                )
            )

            # Rotor reaction: -F_c on the groove at O_g, lever (O_g - O_b) = l_g * n.
            curved_reaction_x -= fc_x
            curved_reaction_y -= fc_y
            rotor_curved_moment += (lever * n_hat[0]) * (-fc_y) - (lever * n_hat[1]) * (-fc_x)

            if want_films:
                films.append((arc.min_film_thickness_m, slider.min_film_thickness_m))

        # Rotor lateral EOM: m_r e_j'' = F_gas + F_journal + sum(-F_c) + m_r omega^2 O_j.
        acc_e_jx = (gx + f_jx + curved_reaction_x + m_r * omega * omega * pin[0]) / m_r
        acc_e_jy = (gy + f_jy + curved_reaction_y + m_r * omega * omega * pin[1]) / m_r

        # Rotor attitude EOM about O_b: I_r phi_orient'' = M_gas + M_fric + curved moment,
        # with phi_orient'' = phi_ddot_ref + dphi_r'' -> subtract the reference inertia.
        eps_j = min(hypot(e_jx, e_jy) / journal_clearance_m, _MAX_ECCENTRICITY_RATIO)
        omega_rel = journal_relative_speed_rad_s(geometry, theta)
        m_fric = -copysign(eccentric_friction_torque_nm(eps_j, omega_rel), phi_dot)
        acc_dphi_r = (gm + m_fric + rotor_curved_moment) / i_r - phi_ddot_ref

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
            return accels, films
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
    state0 = [
        0.0,  # e_jx
        0.0,  # e_jy
        0.0,  # dphi_r
        shift,  # x_IN = +piece_shift (world)
        groove_y0,  # y_IN
        0.0,  # phi_IN
        -shift,  # x_OUT = -piece_shift
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
    mu = LUBRICANT_VISCOSITY_PA_S
    arc_length = 2.0 * radius * half_arc
    for j in range(samples):
        theta = float(omega * t_eval[j])
        state = y[:, j]
        _, films = evaluate(theta, state, want_films=True)
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

        # Couette bush-film friction (comparable to bush_film.friction_power_cycle_w):
        # curved sliding = (rotor swing - piece spin) * r_b; flat sliding = piece world-y.
        flat_length = flat_contact_length_m(geometry, theta, bush)
        total = 0.0
        for ip, iy, cfilm, ffilm in (
            (5, 4, in_curved[j], in_flat[j]),
            (8, 7, out_curved[j], out_flat[j]),
        ):
            u_curved = abs(phi_dot - state[ip + 9]) * radius
            total += mu * u_curved * u_curved / cfilm * arc_length * height
            u_flat = abs(state[iy + 9])
            total += mu * u_flat * u_flat / ffilm * flat_length * height
        friction[j] = total

    bush_friction = float(np.mean(friction))
    quasi_static = friction_power_cycle_w(geometry, trace=trace)

    peak_rotor_ecc = float(np.max(np.hypot(y[0], y[1])) / journal_clearance_m)

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
        rotor_mass_kg=m_r,
        rotor_inertia_kg_m2=i_r,
        piece_mass_kg=m_p,
        piece_inertia_kg_m2=i_p,
        journal_clearance_m=journal_clearance_m,
        revolutions=revolutions,
        frequency_hz=geometry.frequency_hz,
    )
