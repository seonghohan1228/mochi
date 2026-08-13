"""Crank-pin journal bearing friction (Petroff hydrodynamic model).

The crank eccentric journal carries the rotor reaction ``R_j`` (Section 4.6, peak
~3.3 kN) and is the dominant mechanical loss. This adds its friction torque to
the drive torque using Petroff's concentric full-film law with the same POE VG68
oil as the swing-bush film (Section 3.6). See PHYSICS.md Section 4.7.

The oil is sheared at the **relative** angular speed between the crank pin
(rotating at the shaft speed) and the rotor bore. The rotor only swings
``+/-10.37`` deg about 90 deg, so its spin ``omega_rotor = omega * dphi/dtheta``
is a small oscillation and ``omega_rel = omega * (1 - dphi/dtheta)`` stays near
the shaft speed (~154-222 rad/s).

Petroff's law (concentric journal, full film, load-independent):

    T_j = 2 pi mu |omega_rel| r_j**3 L_j / c_j,

with ``r_j`` the journal radius, ``L_j`` its length, ``c_j`` the radial
clearance, ``mu`` the oil viscosity. The dissipation is ``P_j = T_j |omega_rel|``;
referred to the shaft it adds ``P_j / omega`` to the drive torque. The Sommerfeld
number is reported as a validity indicator (Petroff is the concentric
``S -> infinity`` limit; the torque itself does not depend on the load).

The journal radius and length are taken from the CAD model (``r_j = 14.2 mm``,
``L_j = 21 mm``); the radial clearance is **not** drawn (the CAD shows a nominal
full-contact fit), so it stays an assumed ``15 micron`` oil-film value, exposed
for the caller to override. The radius clears the swing-bush groove,
``r_j < cutout_offset - cutout_radius = 17 mm``.
"""

from dataclasses import dataclass
from math import inf, isfinite, pi

from mochi.bush_film import LUBRICANT_VISCOSITY_PA_S
from mochi.kinematics import RotaryGeometry, prescribed_state

# Crank-pin journal geometry. Radius and length are from the CAD model; they
# clear the swing-bush groove (cutout_offset - cutout_radius = 25 - 8 = 17 mm).
# The radial clearance is not drawn (nominal full-contact fit), so it stays an
# assumed oil-film value; override as needed.
JOURNAL_RADIUS_M = 14.2e-3
JOURNAL_LENGTH_M = 21.0e-3
JOURNAL_CLEARANCE_M = 15.0e-6

# Crank-angle step for the central-difference rotor spin rate.
_SPIN_STEP_RAD = 1.0e-5


@dataclass(frozen=True, slots=True)
class JournalFriction:
    """Crank-pin journal friction at one crank angle (SI units)."""

    crank_angle_rad: float
    relative_speed_rad_s: float
    friction_torque_nm: float
    friction_power_w: float
    sommerfeld_number: float


def journal_relative_speed_rad_s(geometry: RotaryGeometry, crank_angle_rad: float) -> float:
    """Relative sliding speed at the crank-pin/rotor journal, ``omega*(1 - dphi/dtheta)``.

    The crank pin turns at the shaft speed; the rotor bore turns at the rotor
    spin ``omega*dphi/dtheta``, so the oil shears at their difference. ``dphi/
    dtheta`` is the central difference of the rotor orientation (the +/-10.37 deg
    swing), matching the sliding-velocity primitive of :mod:`mochi.bush_film`.
    """

    def orientation(angle: float) -> float:
        return prescribed_state(geometry, angle).rotor_orientation_rad

    d_orientation = (
        orientation(crank_angle_rad + _SPIN_STEP_RAD)
        - orientation(crank_angle_rad - _SPIN_STEP_RAD)
    ) / (2.0 * _SPIN_STEP_RAD)
    return geometry.angular_speed_rad_s * (1.0 - d_orientation)


def petroff_friction(
    geometry: RotaryGeometry,
    crank_angle_rad: float,
    *,
    load_n: float | None = None,
    viscosity_pa_s: float = LUBRICANT_VISCOSITY_PA_S,
    radius_m: float = JOURNAL_RADIUS_M,
    length_m: float = JOURNAL_LENGTH_M,
    clearance_m: float = JOURNAL_CLEARANCE_M,
) -> JournalFriction:
    """Return the Petroff journal friction torque, power, and Sommerfeld number.

    ``T_j = 2 pi mu |omega_rel| r**3 L / c`` (concentric, load-independent) and
    ``P_j = T_j |omega_rel|``. If ``load_n`` (the journal reaction magnitude
    ``|R_j|``) is given, the Sommerfeld number ``S = (mu N / P)(r/c)**2`` is
    returned with ``N = |omega_rel| / (2 pi)`` and projected pressure
    ``P = load / (2 r L)``; otherwise it is ``inf`` (no load, concentric limit).
    """

    if not isfinite(crank_angle_rad):
        raise ValueError("Crank angle must be finite.")
    if not isfinite(viscosity_pa_s) or viscosity_pa_s <= 0.0:
        raise ValueError("Viscosity must be a positive, finite value in Pa*s.")
    for name, value in (("radius", radius_m), ("length", length_m), ("clearance", clearance_m)):
        if not isfinite(value) or value <= 0.0:
            raise ValueError(f"Journal {name} must be a positive, finite length in metres.")
    geometry.validate()

    relative_speed = journal_relative_speed_rad_s(geometry, crank_angle_rad)
    speed = abs(relative_speed)
    torque = 2.0 * pi * viscosity_pa_s * speed * radius_m**3 * length_m / clearance_m
    power = torque * speed

    sommerfeld = inf
    if load_n is not None:
        if not isfinite(load_n) or load_n < 0.0:
            raise ValueError("Journal load must be a non-negative, finite force in newtons.")
        if load_n > 0.0:
            revolutions_per_s = speed / (2.0 * pi)
            projected_pressure = load_n / (2.0 * radius_m * length_m)
            sommerfeld = (
                viscosity_pa_s
                * revolutions_per_s
                / projected_pressure
                * (radius_m / clearance_m) ** 2
            )

    return JournalFriction(
        crank_angle_rad=crank_angle_rad,
        relative_speed_rad_s=relative_speed,
        friction_torque_nm=torque,
        friction_power_w=power,
        sommerfeld_number=sommerfeld,
    )
