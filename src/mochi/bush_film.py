"""Swing-bush lubrication films as incompressible Couette + Poiseuille flow.

Each swing-bush piece runs on two ~10 micrometre oil films: a **curved** film
between the bush outside diameter and the rotor's circular groove, and a
**flat** film between the bush flat and the fixed vane. The bush holds a fixed
attitude and only translates with the groove centre (PHYSICS.md section 3.4),
so the two films are driven by different no-slip wall motions:

* flat film: the bush translates along the fixed vane -> sliding speed is the
  groove-centre translation rate (up to ~0.86 m/s at 30 Hz);
* curved film: the rotor tilts +/-10.37 deg relative to the fixed-attitude bush
  -> sliding speed is that orientation-swing rate times the bush radius
  (up to ~0.27 m/s), **not** the full shaft speed.

Each film is treated as a 1-D incompressible thin film of uniform gap ``h``:
Couette shear from the wall motion plus a Poiseuille pressure that is set by the
gas pressures at the two film ends (the chamber pressure at the mouth/OD end and
the 4.0 MPa recess pressure at the inner end). With a uniform gap the pressure
is linear between the two end pressures; the wall shear traction is ``mu U / h``.
The film also reported: leakage flux (Couette + Poiseuille) and friction power.

Lubricant: POE ISO VG68 at ~80 C, ``mu = 0.010 Pa*s`` (neat oil, refrigerant
dilution neglected). See PHYSICS.md section 3.6. This is a quasi-static film
model on the prescribed motion; it carries no hydrodynamic (wedge) load, which
would need a Reynolds/eccentricity model.
"""

from dataclasses import dataclass
from math import isfinite, pi, sqrt

from mochi.chamber_volume import SwingBush
from mochi.chambers import RECESS_PRESSURE_PA, CycleTrace, port_timed_pressures
from mochi.kinematics import RotaryGeometry, prescribed_state, vane_tip_round

# POE ISO VG68 dynamic viscosity at ~80 C operating temperature (Pa*s).
LUBRICANT_VISCOSITY_PA_S = 0.010

# Axial running clearance per side, so the bush height is cylinder height minus
# both; sets the film "width" (the dimension along the bush axis).
AXIAL_CLEARANCE_M = 0.0085e-3

# Crank-angle step for the central-difference sliding velocities.
_VELOCITY_STEP_RAD = 1.0e-5


@dataclass(frozen=True, slots=True)
class FilmFace:
    """One lubrication film (flat or curved) on one bush piece (SI units).

    ``inlet_pressure_pa`` is the chamber-end gas pressure and
    ``outlet_pressure_pa`` the recess-end pressure; with a uniform gap the film
    pressure runs linearly between them. Positive ``slide_velocity`` and
    ``shear_stress`` share the same sign; ``friction_power_w`` is dissipative
    (non-negative). ``leakage_m3_s`` is signed toward increasing film coordinate.
    """

    name: str
    length_m: float
    film_thickness_m: float
    slide_velocity_m_s: float
    shear_stress_pa: float
    inlet_pressure_pa: float
    outlet_pressure_pa: float
    normal_force_n: float
    friction_force_n: float
    leakage_m3_s: float
    friction_power_w: float

    def pressure_at(self, fraction: float) -> float:
        """Linear film pressure at ``fraction`` in ``[0, 1]`` along the film."""

        return (
            self.inlet_pressure_pa + (self.outlet_pressure_pa - self.inlet_pressure_pa) * fraction
        )


@dataclass(frozen=True, slots=True)
class BushPieceFilm:
    """The flat and curved films of one bush piece (the IN or OUT side)."""

    side: str
    flat: FilmFace
    curved: FilmFace


@dataclass(frozen=True, slots=True)
class BushFilmState:
    """Both bush pieces' films at one crank angle, with bush totals."""

    crank_angle_rad: float
    viscosity_pa_s: float
    in_piece: BushPieceFilm
    out_piece: BushPieceFilm
    total_friction_power_w: float
    total_leakage_m3_s: float


def _central_difference(function, angle_rad: float, step_rad: float) -> float:
    return (function(angle_rad + step_rad) - function(angle_rad - step_rad)) / (2.0 * step_rad)


def flat_slide_velocity(geometry: RotaryGeometry, crank_angle_rad: float) -> float:
    """Bush translation speed along the vane (flat-film driver), m/s.

    The bush centre follows the groove centre ``cutout_center_m``; its vane-axis
    translation rate is ``omega * d(cutout_center_y)/d(theta)``.
    """

    def cutout_y(angle: float) -> float:
        return prescribed_state(geometry, angle).cutout_center_m[1]

    d_y = _central_difference(cutout_y, crank_angle_rad, _VELOCITY_STEP_RAD)
    return geometry.angular_speed_rad_s * d_y


def curved_slide_velocity(geometry: RotaryGeometry, crank_angle_rad: float) -> float:
    """Rotor-vs-bush sliding speed on the curved film, m/s.

    The fixed-attitude bush translates with the groove centre, so only the
    rotor's orientation swing shears the curved film:
    ``omega * d(rotor_orientation)/d(theta) * r_bush``.
    """

    def orientation(angle: float) -> float:
        return prescribed_state(geometry, angle).rotor_orientation_rad

    d_orientation = _central_difference(orientation, crank_angle_rad, _VELOCITY_STEP_RAD)
    return geometry.angular_speed_rad_s * d_orientation * SwingBush().piece_outer_radius_m


def _film_face(
    name: str,
    length_m: float,
    thickness_m: float,
    slide_velocity_m_s: float,
    height_m: float,
    viscosity_pa_s: float,
    inlet_pressure_pa: float,
    outlet_pressure_pa: float,
) -> FilmFace:
    shear_stress = viscosity_pa_s * slide_velocity_m_s / thickness_m
    pressure_gradient = (outlet_pressure_pa - inlet_pressure_pa) / length_m
    # Flux per unit width: Couette drag minus Poiseuille back-flow.
    flux_per_width = (
        0.5 * slide_velocity_m_s * thickness_m
        - (thickness_m**3 / (12.0 * viscosity_pa_s)) * pressure_gradient
    )
    return FilmFace(
        name=name,
        length_m=length_m,
        film_thickness_m=thickness_m,
        slide_velocity_m_s=slide_velocity_m_s,
        shear_stress_pa=shear_stress,
        inlet_pressure_pa=inlet_pressure_pa,
        outlet_pressure_pa=outlet_pressure_pa,
        normal_force_n=0.5 * (inlet_pressure_pa + outlet_pressure_pa) * length_m * height_m,
        friction_force_n=shear_stress * length_m * height_m,
        leakage_m3_s=flux_per_width * height_m,
        friction_power_w=abs(shear_stress * slide_velocity_m_s) * length_m * height_m,
    )


def film_thicknesses_m(geometry: RotaryGeometry, bush: SwingBush) -> tuple[float, float]:
    """Return ``(flat_film, curved_film)`` gaps from the shift geometry (m).

    Both are 0.010 mm for the supplied geometry: the flat film is
    ``flat_offset + shift - vane_half`` and the curved film is
    ``cutout_radius - bush_radius - shift`` (PHYSICS.md section 3.3).
    """

    flat_film = bush.flat_offset_m + bush.piece_shift_m - 0.5 * geometry.vane_width_m
    curved_film = geometry.cutout_radius_m - bush.piece_outer_radius_m - bush.piece_shift_m
    return flat_film, curved_film


def full_flat_contact_length_m(bush: SwingBush) -> float:
    """Full bush flat-face contact chord length against the vane flank (m)."""

    return 2.0 * sqrt(
        (bush.piece_outer_radius_m - bush.fillet_radius_m) ** 2
        - (bush.flat_offset_m + bush.fillet_radius_m) ** 2
    )


def flat_contact_length_m(
    geometry: RotaryGeometry,
    crank_angle_rad: float,
    bush: SwingBush | None = None,
) -> float:
    """Effective bush flat-face contact length on the vane flank at one crank angle.

    The bush flat spans the groove centre plus/minus half the full chord, but only
    the **straight** vane flank above the tip round carries it (``y >= y_tip + R``,
    :func:`mochi.kinematics.vane_tip_round`). Near BDC the groove centre drops so
    the lowest part of the flat overhangs the R1.5 tip corner and the contact
    shortens from the full 11.94 mm to about 11.47 mm; it is full away from BDC.
    """

    bush = SwingBush() if bush is None else bush
    full_flat = full_flat_contact_length_m(bush)
    half_flat = 0.5 * full_flat
    groove_y = prescribed_state(geometry, crank_angle_rad).cutout_center_m[1]
    _, flank_top, _, _ = vane_tip_round(geometry)
    lower = max(groove_y - half_flat, flank_top)
    return max(min(groove_y + half_flat - lower, full_flat), 0.0)


def film_state(
    geometry: RotaryGeometry,
    crank_angle_rad: float,
    *,
    viscosity_pa_s: float = LUBRICANT_VISCOSITY_PA_S,
    bush: SwingBush | None = None,
    trace: CycleTrace | None = None,
    in_chamber_pressure_pa: float | None = None,
    out_chamber_pressure_pa: float | None = None,
    recess_pressure_pa: float = RECESS_PRESSURE_PA,
) -> BushFilmState:
    """Solve both bush pieces' films at one crank angle.

    The IN piece faces the suction chamber and the OUT piece the compression
    chamber; each film runs between that chamber pressure (mouth/OD end) and the
    recess pressure (inner end). Chamber pressures default to the port-timed
    rule (section 3.4) but can be overridden. Boundary-condition assignment is a
    first-cut modelling choice (the exact film-edge-to-gas-space topology is a
    deferred open item), hence the explicit pressure arguments.
    """

    if not isfinite(crank_angle_rad):
        raise ValueError("Crank angle must be finite.")
    if not isfinite(viscosity_pa_s) or viscosity_pa_s <= 0.0:
        raise ValueError("Viscosity must be a positive, finite value in Pa*s.")
    geometry.validate()
    bush = SwingBush() if bush is None else bush
    bush.validate(geometry)

    flat_h, curved_h = film_thicknesses_m(geometry, bush)
    if flat_h <= 0.0 or curved_h <= 0.0:
        raise ValueError("Film thicknesses must be positive; check the bush shift geometry.")

    # The flat contact shortens near BDC as the bush overhangs the R1.5 vane-tip
    # round; the curved face (on the rotor groove) is unaffected.
    length_flat = flat_contact_length_m(geometry, crank_angle_rad, bush)
    length_curved = bush.piece_outer_radius_m * 2.0 * bush.half_arc_rad()
    height = geometry.cylinder_height_m - 2.0 * AXIAL_CLEARANCE_M

    velocity_flat = flat_slide_velocity(geometry, crank_angle_rad)
    velocity_curved = curved_slide_velocity(geometry, crank_angle_rad)

    if in_chamber_pressure_pa is None or out_chamber_pressure_pa is None:
        pressures = port_timed_pressures(geometry, crank_angle_rad, trace=trace)
        if in_chamber_pressure_pa is None:
            in_chamber_pressure_pa = pressures.suction_pressure_pa
        if out_chamber_pressure_pa is None:
            out_chamber_pressure_pa = pressures.compression_pressure_pa

    def piece(side: str, chamber_pressure_pa: float) -> BushPieceFilm:
        return BushPieceFilm(
            side=side,
            flat=_film_face(
                "flat",
                length_flat,
                flat_h,
                velocity_flat,
                height,
                viscosity_pa_s,
                chamber_pressure_pa,
                recess_pressure_pa,
            ),
            curved=_film_face(
                "curved",
                length_curved,
                curved_h,
                velocity_curved,
                height,
                viscosity_pa_s,
                chamber_pressure_pa,
                recess_pressure_pa,
            ),
        )

    in_piece = piece("IN", in_chamber_pressure_pa)
    out_piece = piece("OUT", out_chamber_pressure_pa)

    faces = (in_piece.flat, in_piece.curved, out_piece.flat, out_piece.curved)
    return BushFilmState(
        crank_angle_rad=crank_angle_rad,
        viscosity_pa_s=viscosity_pa_s,
        in_piece=in_piece,
        out_piece=out_piece,
        total_friction_power_w=sum(face.friction_power_w for face in faces),
        total_leakage_m3_s=sum(abs(face.leakage_m3_s) for face in faces),
    )


def friction_power_cycle_w(
    geometry: RotaryGeometry,
    *,
    samples: int = 360,
    viscosity_pa_s: float = LUBRICANT_VISCOSITY_PA_S,
    trace: CycleTrace | None = None,
) -> float:
    """Crank-angle-averaged total bush-film friction power over one revolution.

    A mechanical loss term to later subtract from the indicated power for a
    shaft-power estimate.
    """

    if samples < 8:
        raise ValueError("The friction-power average needs at least eight samples.")

    total = 0.0
    for index in range(samples):
        angle = 2.0 * pi * index / samples
        total += film_state(
            geometry, angle, viscosity_pa_s=viscosity_pa_s, trace=trace
        ).total_friction_power_w
    return total / samples
