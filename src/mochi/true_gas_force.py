"""True-geometry gas force and shaft torque (mouth-aware, consistent with the volume).

Section 4.5's :mod:`mochi.gas_force` integrates the gas pressure over the rotor
*outside-diameter disc* -- a circular arc whose closed form makes the resultant
pass exactly through the rotor centre ``O_r`` and the torque ``T_gas = O_r x F``.
That circular-rotor basis is ~3% below the true indicated work (715 W vs 738 W)
because it discards the rotor **mouth** (the recess that admits the vane and
holds the swing bushes; :mod:`mochi.rotor_profile`), the very cavity the volume
model of :mod:`mochi.chamber_volume` keeps.

This module integrates the gas load over the **real rotor contour** instead, so
the force and the torque use the same shape the pressure and volume already use:

* **Force** -- the pressure integral over the true material boundary (OD arc plus
  the mouth walls), skipping spans sealed by a bush piece or covered by the vane.
  A grid-converged raster cross-check confirms it to <1%; the mouth cavity lowers
  the net rotor load ~20% below the circular closed form near peak.

* **Torque** -- on the true geometry the resultant no longer passes through
  ``O_r`` and the rotor both orbits and swings, so ``O_r x F`` is *not* the shaft
  torque, and the true indicated work includes work done on the moving bush too.
  The work-conjugate (shaft) gas torque is therefore taken from virtual work,
  ``T_gas(theta) = -(p_c dV_c/dtheta + p_s dV_s/dtheta)`` on the true-geometry
  volume trace, whose crank-angle integral equals the true indicated work
  (738 W) exactly -- the cross-check that binds this rung to Section 3.5.

The circular closed form of Section 4.5 is kept as the fast analytic check; this
module is the refined headline for the gas force, the crank-pin reaction
``R_j = -F_gas`` (Section 4.6), and the indicated gas torque (Section 4.5).
"""

from dataclasses import dataclass
from math import degrees, hypot, isfinite, pi

import numpy as np

from mochi.chamber_volume import SwingBush
from mochi.chambers import (
    CycleTrace,
    build_cycle_trace,
    in_seal_over_window,
    port_timed_pressures,
    seal_over_half_angle_rad,
)
from mochi.indicated_work import indicated_work_j
from mochi.kinematics import RotaryGeometry, prescribed_state
from mochi.rotor_profile import MouthGeometry, rotor_contour

FULL_TURN_RAD = 2.0 * pi
# Central-difference half-step for dV/dtheta; small against the trace spacing.
VOLUME_DERIVATIVE_STEP_RAD = 1.0e-4


@dataclass(frozen=True, slots=True)
class TrueGasLoad:
    """Mouth-aware gas load at one crank angle (global frame, SI units)."""

    crank_angle_rad: float
    rotor_force_n: tuple[float, float]
    rotor_force_mag_n: float
    rotor_torque_nm: float
    rotor_moment_about_center_nm: float
    suction_pressure_pa: float
    compression_pressure_pa: float


@dataclass(frozen=True, slots=True)
class TrueGasTorqueWork:
    """True-geometry gas-torque work cross-check (SI units)."""

    torque_work_j: float
    power_w: float
    indicated_work_j: float
    indicated_power_w: float
    circular_power_w: float
    frequency_hz: float
    samples: int


def _wetted_rotor_force(
    geometry: RotaryGeometry,
    crank_angle_rad: float,
    suction_pa: float,
    compression_pa: float,
    mouth: MouthGeometry | None,
    bush: SwingBush,
) -> tuple[float, float, float]:
    """Net gas force and moment on the true rotor boundary at one crank angle.

    Integrates ``dF = -p n ds H`` over the closed rotor material contour, with
    the outward normal from the traversal orientation. A segment is skipped where
    its midpoint is covered by a bush piece or lies under the vane column; the
    remaining spans take the suction or compression pressure by the bore angle of
    their midpoint relative to the rotor-cylinder contact (bore angle ``theta``).
    Returns ``(F_x, F_y, M_z)`` with ``M_z`` the moment about the rotor centre
    ``O_b`` (the axis the rotor attitude EOM of the D2 dynamics uses).
    """

    state = prescribed_state(geometry, crank_angle_rad)
    groove = state.cutout_center_m
    vane_tip_y = state.vane_tip_m[1]
    half_width = 0.5 * geometry.vane_width_m
    height = geometry.cylinder_height_m
    contact_deg = degrees(crank_angle_rad) % 360.0

    loop = np.asarray(rotor_contour(geometry, crank_angle_rad, mouth=mouth).material, dtype=float)
    start = loop
    end = np.roll(loop, -1, axis=0)
    tangent = end - start
    length = np.hypot(tangent[:, 0], tangent[:, 1])
    good = length > 0.0

    signed_area2 = float(np.sum(start[:, 0] * end[:, 1] - end[:, 0] * start[:, 1]))
    if signed_area2 > 0.0:  # CCW: outward normal is the tangent turned -90 deg
        normal = np.stack([tangent[:, 1], -tangent[:, 0]], axis=1)
    else:
        normal = np.stack([-tangent[:, 1], tangent[:, 0]], axis=1)
    normal = np.divide(normal, np.where(good, length, 1.0)[:, None])

    mid = 0.5 * (start + end)
    covered = bush.occupies(mid[:, 0], mid[:, 1], groove) | (
        (np.abs(mid[:, 0]) <= half_width) & (mid[:, 1] >= vane_tip_y)
    )
    bore_deg = np.degrees(np.arctan2(mid[:, 0], mid[:, 1])) % 360.0
    pressure = np.where(bore_deg < contact_deg, suction_pa, compression_pa)
    pressure = np.where(covered | ~good, 0.0, pressure)

    scale = -pressure * length * height
    d_force_x = scale * normal[:, 0]
    d_force_y = scale * normal[:, 1]
    force_x = float(np.sum(d_force_x))
    force_y = float(np.sum(d_force_y))
    # Moment about the rotor centre O_b: sum (r - O_b) x dF over the contour.
    rotor_x, rotor_y = state.rotor_center_m
    arm_x = mid[:, 0] - rotor_x
    arm_y = mid[:, 1] - rotor_y
    moment_z = float(np.sum(arm_x * d_force_y - arm_y * d_force_x))
    return force_x, force_y, moment_z


def _virtual_torque_nm(
    geometry: RotaryGeometry,
    crank_angle_rad: float,
    trace: CycleTrace,
    suction_pa: float,
    compression_pa: float,
    step_rad: float = VOLUME_DERIVATIVE_STEP_RAD,
) -> float:
    """Work-conjugate shaft gas torque ``-(p_c dV_c/dtheta + p_s dV_s/dtheta)``."""

    d_compression = (
        trace.compression_at(crank_angle_rad + step_rad)
        - trace.compression_at(crank_angle_rad - step_rad)
    ) / (2.0 * step_rad)
    d_suction = (
        trace.suction_at(crank_angle_rad + step_rad) - trace.suction_at(crank_angle_rad - step_rad)
    ) / (2.0 * step_rad)
    return -(compression_pa * d_compression + suction_pa * d_suction)


def true_gas_load(
    geometry: RotaryGeometry,
    crank_angle_rad: float,
    *,
    trace: CycleTrace | None = None,
    mouth: MouthGeometry | None = None,
    bush: SwingBush | None = None,
) -> TrueGasLoad:
    """Return the mouth-aware gas force and shaft gas torque at one crank angle.

    The force is the pressure integral over the true rotor contour (mouth-aware);
    the torque is the virtual-work shaft torque on the true-geometry volume trace.
    Raises ``ValueError`` inside the seal-over window (no separated contact) and
    for a non-finite crank angle.
    """

    geometry.validate()
    if not isfinite(crank_angle_rad):
        raise ValueError("Crank angle must be finite.")
    if in_seal_over_window(geometry, crank_angle_rad):
        raise ValueError(
            "True gas load is undefined in the seal-over window: the rotor-cylinder "
            "contact hides under the vane and the chambers merge."
        )

    trace = build_cycle_trace(geometry) if trace is None else trace
    bush = SwingBush() if bush is None else bush
    pressures = port_timed_pressures(geometry, crank_angle_rad, trace=trace)
    suction_pa = pressures.suction_pressure_pa
    compression_pa = pressures.compression_pressure_pa

    fx, fy, moment = _wetted_rotor_force(
        geometry, crank_angle_rad, suction_pa, compression_pa, mouth, bush
    )
    torque = _virtual_torque_nm(geometry, crank_angle_rad, trace, suction_pa, compression_pa)

    return TrueGasLoad(
        crank_angle_rad=crank_angle_rad,
        rotor_force_n=(fx, fy),
        rotor_force_mag_n=hypot(fx, fy),
        rotor_torque_nm=torque,
        rotor_moment_about_center_nm=moment,
        suction_pressure_pa=suction_pa,
        compression_pressure_pa=compression_pa,
    )


def true_gas_torque_work_j(
    geometry: RotaryGeometry,
    *,
    samples: int = 720,
    trace: CycleTrace | None = None,
) -> TrueGasTorqueWork:
    """Cross-check the true shaft gas torque against the indicated work.

    Integrates the virtual-work torque over the same range as
    :func:`mochi.indicated_work.indicated_work_j` (``[0, 2*pi - alpha]``, the
    seal-over entry), so ``torque_work_j`` equals the true indicated work exactly
    and ``power_w`` reproduces the 738 W headline -- above the circular-rotor
    715 W of Section 4.5.
    """

    geometry.validate()
    if samples < 8:
        raise ValueError("The true gas-torque cross-check needs at least eight samples.")

    trace = build_cycle_trace(geometry) if trace is None else trace
    seal_over_entry = FULL_TURN_RAD - seal_over_half_angle_rad(geometry)
    angles = tuple(seal_over_entry * index / samples for index in range(samples + 1))
    pressures = tuple(port_timed_pressures(geometry, angle, trace=trace) for angle in angles)
    compression_p = tuple(p.compression_pressure_pa for p in pressures)
    suction_p = tuple(p.suction_pressure_pa for p in pressures)
    compression_v = tuple(trace.compression_at(angle) for angle in angles)
    suction_v = tuple(trace.suction_at(angle) for angle in angles)

    # oint T dtheta = -oint (p_c dV_c + p_s dV_s): evaluated as -sum p dV so it
    # telescopes cleanly across the 0/2*pi seam (differentiating the trace there
    # would pick up the seal-over volume jump). This is exactly the identity that
    # binds the true shaft torque to the indicated work of Section 3.5.
    torque_work = 0.0
    for index in range(samples):
        torque_work -= (
            0.5
            * (compression_p[index] + compression_p[index + 1])
            * (compression_v[index + 1] - compression_v[index])
        )
        torque_work -= (
            0.5
            * (suction_p[index] + suction_p[index + 1])
            * (suction_v[index + 1] - suction_v[index])
        )

    frequency = geometry.frequency_hz
    indicated = indicated_work_j(geometry, trace=trace, samples=samples)
    # The circular-rotor closed-form basis (Section 4.5), for the loss-ladder note.
    from mochi.gas_force import gas_torque_work_j

    circular = gas_torque_work_j(geometry, trace=trace)
    return TrueGasTorqueWork(
        torque_work_j=torque_work,
        power_w=torque_work * frequency,
        indicated_work_j=indicated.net_work_j,
        indicated_power_w=indicated.power_w,
        circular_power_w=circular.power_w,
        frequency_hz=frequency,
        samples=samples,
    )


def peak_rotor_force_n(
    geometry: RotaryGeometry,
    *,
    samples: int = 360,
    trace: CycleTrace | None = None,
    mouth: MouthGeometry | None = None,
    bush: SwingBush | None = None,
) -> float:
    """Return the peak true-geometry rotor gas load over one revolution.

    Swept over the separated range ``[alpha, 2*pi - alpha]`` (outside seal-over),
    this is the peak crank-pin journal reaction ``|R_j| = |F_gas|`` (Section 4.6).
    """

    geometry.validate()
    if samples < 8:
        raise ValueError("The peak-force sweep needs at least eight samples.")

    trace = build_cycle_trace(geometry) if trace is None else trace
    bush = SwingBush() if bush is None else bush
    half_angle = seal_over_half_angle_rad(geometry)
    span = FULL_TURN_RAD - 2.0 * half_angle
    peak = 0.0
    for index in range(samples + 1):
        angle = half_angle + span * index / samples
        pressures = port_timed_pressures(geometry, angle, trace=trace)
        fx, fy, _moment = _wetted_rotor_force(
            geometry,
            angle,
            pressures.suction_pressure_pa,
            pressures.compression_pressure_pa,
            mouth,
            bush,
        )
        peak = max(peak, hypot(fx, fy))
    return peak
