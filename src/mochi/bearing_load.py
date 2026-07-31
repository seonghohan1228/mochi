"""Crank-pin bearing reaction and drive torque along the prescribed motion.

This closes the rotor's quasi-static force balance (net DOF = 1, ``theta``
prescribed) on top of the gas force of Section 4.5: the crank eccentric
(journal) bearing reaction that carries the gas load, and the drive torque the
shaft must supply. See PHYSICS.md Section 4.6.

**Statical indeterminacy (resolved by the standard treatment).** The crank
journal fixes the rotor centre on the line ``O`-`O_r``, and the rotor is ideally
tangent to the cylinder at ``C`` with ``O, O_r, C`` colinear, so the
rotor-cylinder contact normal force and the radial journal reaction lie on the
*same* line -- rigid-body statics cannot split them. A moment balance about the
rotor centre ``O_r`` kills both (zero moment arm) and leaves only friction
torques, so no independent contact normal force exists under the rigid,
frictionless, inertialess model. The accepted treatment (Yanagisawa & Shimizu
1982/1985; Aw & Ooi 2021 review) lumps the radial load onto the journal bearing
and treats the contact as friction-only. In the frictionless quasi-static limit
the journal therefore carries the whole gas force,

    R_j = -F_gas,

and the contact normal force is not a determinate rigid-body quantity (it would
require an oil-film / Hertzian compliance model, a later rung).

**Drive torque.** ``T_drive = T_gas + T_friction``: the indicated gas torque of
Section 4.5 plus the mechanical friction referred to the shaft,
``T_friction = (P_bush + P_journal) / omega``. Two friction terms are modelled:
the swing-bush film dissipation (Section 3.6) and the crank-pin journal bearing
friction (Section 4.7, Petroff); the latter dominates. Dry rotor-cylinder
contact friction is still omitted (indeterminate contact force, no coefficient).

Cross-check: the cycle-mean shaft power equals the indicated power (Section 3.5)
plus the bush-film and journal friction powers (Sections 3.6, 4.7), binding this
rung to the gas force (Section 4.5) and the earlier physics rungs.
"""

from dataclasses import dataclass
from math import hypot, isfinite, pi

from mochi.bush_film import LUBRICANT_VISCOSITY_PA_S, film_state
from mochi.chambers import CycleTrace, build_cycle_trace, seal_over_half_angle_rad
from mochi.gas_force import gas_load
from mochi.indicated_work import indicated_work_j
from mochi.journal_bearing import petroff_friction
from mochi.kinematics import RotaryGeometry
from mochi.true_gas_force import peak_rotor_force_n

FULL_TURN_RAD = 2.0 * pi


@dataclass(frozen=True, slots=True)
class MechanismLoad:
    """Crank-pin reaction and drive torque at one crank angle (global frame, SI)."""

    crank_angle_rad: float
    journal_force_n: tuple[float, float]
    journal_force_mag_n: float
    gas_torque_nm: float
    bush_friction_torque_nm: float
    journal_friction_torque_nm: float
    friction_torque_nm: float
    drive_torque_nm: float


@dataclass(frozen=True, slots=True)
class ShaftWork:
    """Drive-torque cross-check against indicated + friction power (SI units)."""

    drive_work_j: float
    gas_work_j: float
    friction_work_j: float
    shaft_power_w: float
    indicated_power_w: float
    friction_power_w: float
    bush_friction_power_w: float
    journal_friction_power_w: float
    peak_journal_force_n: float
    frequency_hz: float
    samples: int


def mechanism_load(
    geometry: RotaryGeometry,
    crank_angle_rad: float,
    *,
    trace: CycleTrace | None = None,
    viscosity_pa_s: float = LUBRICANT_VISCOSITY_PA_S,
) -> MechanismLoad:
    """Return the crank-pin bearing reaction and drive torque at one crank angle.

    The journal reaction is ``R_j = -F_gas`` (the radial gas load lumped onto the
    journal, per Section 4.6); the drive torque is the gas torque plus the
    bush-film (Section 3.6) and journal (Section 4.7) friction referred to the
    shaft. Raises ``ValueError`` inside the seal-over window (via
    :func:`mochi.gas_force.gas_load`) and for a non-finite angle.
    """

    geometry.validate()
    if not isfinite(crank_angle_rad):
        raise ValueError("Crank angle must be finite.")

    # gas_load rejects the seal-over window before touching the volume trace, so
    # pass the caller's trace straight through and let it build only when valid.
    load = gas_load(geometry, crank_angle_rad, trace=trace)
    journal = (-load.rotor_force_n[0], -load.rotor_force_n[1])
    journal_magnitude = hypot(journal[0], journal[1])

    shaft_speed = geometry.angular_speed_rad_s
    bush_power_w = film_state(
        geometry, crank_angle_rad, viscosity_pa_s=viscosity_pa_s, trace=trace
    ).total_friction_power_w
    journal_power_w = petroff_friction(
        geometry, crank_angle_rad, load_n=journal_magnitude, viscosity_pa_s=viscosity_pa_s
    ).friction_power_w
    bush_torque = bush_power_w / shaft_speed
    journal_torque = journal_power_w / shaft_speed
    friction_torque = bush_torque + journal_torque

    return MechanismLoad(
        crank_angle_rad=crank_angle_rad,
        journal_force_n=journal,
        journal_force_mag_n=journal_magnitude,
        gas_torque_nm=load.rotor_torque_nm,
        bush_friction_torque_nm=bush_torque,
        journal_friction_torque_nm=journal_torque,
        friction_torque_nm=friction_torque,
        drive_torque_nm=load.rotor_torque_nm + friction_torque,
    )


def shaft_work_j(
    geometry: RotaryGeometry,
    *,
    samples: int = 720,
    trace: CycleTrace | None = None,
    viscosity_pa_s: float = LUBRICANT_VISCOSITY_PA_S,
) -> ShaftWork:
    """Integrate the drive torque and return the shaft power over a revolution.

    ``shaft_power_w = indicated_power_w + friction_power_w`` — the isentropic-ideal
    **baseline** indicated power (Section 3.5) plus the bush (Section 3.6) and
    journal (Section 4.7) friction — so it exceeds the indicated power by the
    friction. ``friction_power_w`` splits into ``bush_friction_power_w`` and
    ``journal_friction_power_w``. The discharge reed-valve overpressure (Section
    3.8, ~51 W) is a **separate performance term not propagated into the loads**;
    it is not added here.

    ``drive_work_j`` integrates ``oint T_drive d(theta)`` and splits exactly into
    ``gas_work_j + friction_work_j`` (a torque-integral cross-check); note
    ``drive_work_j * f`` is **not** the shaft power, because ``gas_work_j`` is the
    gas torque on the circular-rotor basis (Section 4.5), ~3% below the true
    indicated work, and that basis gap would otherwise swamp the friction term and
    make the shaft power spuriously fall below the indicated power.

    Integration runs over the separated range ``[alpha, 2*pi - alpha]`` (the
    seal-over half-angle on each end), where two distinct chambers exist, exactly
    as in :func:`mochi.gas_force.gas_torque_work_j`.
    """

    geometry.validate()
    if samples < 8:
        raise ValueError("The shaft-power cross-check needs at least eight samples.")

    trace = build_cycle_trace(geometry) if trace is None else trace
    half_angle = seal_over_half_angle_rad(geometry)
    start = half_angle
    span = FULL_TURN_RAD - 2.0 * half_angle
    angles = tuple(start + span * index / samples for index in range(samples + 1))
    loads = tuple(
        mechanism_load(geometry, angle, trace=trace, viscosity_pa_s=viscosity_pa_s)
        for angle in angles
    )

    drive_work = 0.0
    gas_work = 0.0
    friction_work = 0.0
    bush_work = 0.0
    journal_work = 0.0
    for index in range(samples):
        step = angles[index + 1] - angles[index]
        drive_work += 0.5 * (loads[index].drive_torque_nm + loads[index + 1].drive_torque_nm) * step
        gas_work += 0.5 * (loads[index].gas_torque_nm + loads[index + 1].gas_torque_nm) * step
        friction_work += (
            0.5 * (loads[index].friction_torque_nm + loads[index + 1].friction_torque_nm) * step
        )
        bush_work += (
            0.5
            * (loads[index].bush_friction_torque_nm + loads[index + 1].bush_friction_torque_nm)
            * step
        )
        journal_work += (
            0.5
            * (
                loads[index].journal_friction_torque_nm
                + loads[index + 1].journal_friction_torque_nm
            )
            * step
        )

    frequency = geometry.frequency_hz
    indicated_power = indicated_work_j(geometry, trace=trace).power_w
    friction_power = friction_work * frequency
    # Shaft power is the isentropic-ideal BASELINE indicated power (§3.5) plus the
    # friction, so it exceeds the indicated power by the friction alone. The reed-
    # valve overpressure (§3.8) is a separate performance term, not propagated
    # here. It is NOT drive_work * f either: the drive-torque work uses the gas
    # torque on the circular-rotor basis (§4.5), ~3% below the true indicated work,
    # and that basis gap would otherwise swamp the smaller friction term.
    return ShaftWork(
        drive_work_j=drive_work,
        gas_work_j=gas_work,
        friction_work_j=friction_work,
        shaft_power_w=indicated_power + friction_power,
        indicated_power_w=indicated_power,
        friction_power_w=friction_power,
        bush_friction_power_w=bush_work * frequency,
        journal_friction_power_w=journal_work * frequency,
        # Peak reaction on the mouth-aware true force (Section 4.5/4.6); the rotor
        # mouth lowers it ~20% below the circular closed form.
        peak_journal_force_n=peak_rotor_force_n(geometry, trace=trace),
        frequency_hz=frequency,
        samples=samples,
    )
