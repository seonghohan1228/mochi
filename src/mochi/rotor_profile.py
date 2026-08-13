"""Rotor boundary contour for the confirmed compressor geometry.

The rotor is not a plain disc: the inlet side of its outside diameter is cut
by a straight, and the mouth that admits the vane is bounded by two
tangent-continuous lips that differ between the inlet and outlet sides. This
module turns those dimensions into drawable polylines; see PHYSICS.md
section 3.3 for the dimensions and their provenance.
"""

from dataclasses import dataclass
from math import atan2, cos, hypot, isfinite, pi, radians, sin, sqrt

from mochi.kinematics import MM, RotaryGeometry, prescribed_state

Point = tuple[float, float]


@dataclass(frozen=True, slots=True)
class MouthGeometry:
    """Mouth and lip dimensions, all angles at the rotor center."""

    flat_start_rad: float
    flat_end_rad: float
    flat_end_radius_m: float
    outlet_tangent_rad: float
    inlet_sweep_rad: float
    outlet_sweep_rad: float
    lip_radius_m: float
    blend_radius_m: float

    @classmethod
    def default(cls) -> "MouthGeometry":
        """Return the confirmed mouth dimensions (2026-07-20)."""

        return cls(
            flat_start_rad=radians(34.8),
            flat_end_rad=radians(13.4),
            flat_end_radius_m=33.657 * MM,
            outlet_tangent_rad=radians(13.0),
            inlet_sweep_rad=radians(102.2),
            outlet_sweep_rad=radians(92.6),
            lip_radius_m=1.5 * MM,
            blend_radius_m=1.0 * MM,
        )

    def validate(self, geometry: RotaryGeometry) -> None:
        """Raise ``ValueError`` when the mouth cannot sit on this rotor."""

        values = (
            self.flat_start_rad,
            self.flat_end_rad,
            self.flat_end_radius_m,
            self.outlet_tangent_rad,
            self.inlet_sweep_rad,
            self.outlet_sweep_rad,
            self.lip_radius_m,
            self.blend_radius_m,
        )
        if not all(isfinite(value) for value in values):
            raise ValueError("All mouth dimensions must be finite.")
        if not 0.0 < self.flat_end_rad < self.flat_start_rad < 0.5 * pi:
            raise ValueError("Inlet flat must run from a larger to a smaller angle.")
        if not 0.0 < self.outlet_tangent_rad < 0.5 * pi:
            raise ValueError("Outlet tangent angle must lie between 0 and 90 degrees.")
        if not 0.0 < self.flat_end_radius_m < geometry.rotor_radius_m:
            raise ValueError("Inlet flat must end inside the rotor outside diameter.")
        if self.lip_radius_m <= 0.0 or self.blend_radius_m <= 0.0:
            raise ValueError("Lip radii must be positive.")
        if self.inlet_sweep_rad <= 0.0 or self.outlet_sweep_rad <= 0.0:
            raise ValueError("Lip sweeps must be positive.")


@dataclass(frozen=True, slots=True)
class RotorContour:
    """Rotor boundary at one crank angle, in global metres.

    The three polylines are the real material edges and share endpoints:
    ``od_arc`` runs from the inlet flat start the long way round to the
    outlet tangent point, ``inlet_flat`` from that same start to the flat's
    end, and ``mouth_path`` from the flat's end through both lips and the
    groove back to the outlet tangent point. ``material`` closes them into a
    single filled outline.
    """

    od_arc: tuple[Point, ...]
    inlet_flat: tuple[Point, ...]
    mouth_path: tuple[Point, ...]
    material: tuple[Point, ...]


def _solve_lip(
    start: Point,
    direction: Point,
    sweep_rad: float,
    offset_m: float,
    groove_radius_m: float,
    lip_radius_m: float,
    blend_radius_m: float,
    arc_points: int = 41,
    blend_points: int = 31,
) -> tuple[list[Point], float]:
    """One lip in rotor-local (transverse, axial) metres.

    Travels from ``start`` along ``direction``: a ``lip_radius_m`` arc of
    ``sweep_rad``, a straight, then a ``blend_radius_m`` arc tangent to the
    groove circle. Returns the points and the groove-blend angle.
    """

    centre = (
        start[0] - lip_radius_m * direction[1],
        start[1] + lip_radius_m * direction[0],
    )
    start_az = atan2(start[1] - centre[1], start[0] - centre[0])
    arc = [
        (
            centre[0] + lip_radius_m * cos(start_az + sweep_rad * i / (arc_points - 1)),
            centre[1] + lip_radius_m * sin(start_az + sweep_rad * i / (arc_points - 1)),
        )
        for i in range(arc_points)
    ]
    heading = start_az + sweep_rad + 0.5 * pi
    along = (cos(heading), sin(heading))
    normal = (-along[1], along[0])
    end = arc[-1]
    to_groove_x = end[0] + blend_radius_m * normal[0]
    to_groove_y = end[1] + blend_radius_m * normal[1] - offset_m
    linear = 2.0 * (to_groove_x * along[0] + to_groove_y * along[1])
    constant = (
        to_groove_x * to_groove_x
        + to_groove_y * to_groove_y
        - (groove_radius_m + blend_radius_m) ** 2
    )
    discriminant = linear * linear - 4.0 * constant
    if discriminant < 0.0:
        raise ValueError("Lip cannot be blended into the circular groove.")
    straight = (-linear - sqrt(discriminant)) / 2.0
    joint = (end[0] + straight * along[0], end[1] + straight * along[1])
    blend_centre = (
        joint[0] + blend_radius_m * normal[0],
        joint[1] + blend_radius_m * normal[1],
    )
    span = hypot(blend_centre[0], blend_centre[1] - offset_m)
    touch = (
        groove_radius_m * blend_centre[0] / span,
        offset_m + groove_radius_m * (blend_centre[1] - offset_m) / span,
    )
    joint_az = atan2(joint[1] - blend_centre[1], joint[0] - blend_centre[0])
    touch_az = atan2(touch[1] - blend_centre[1], touch[0] - blend_centre[0])
    blend_sweep = (touch_az - joint_az) % (2.0 * pi)
    blend = [
        (
            blend_centre[0] + blend_radius_m * cos(joint_az + blend_sweep * i / (blend_points - 1)),
            blend_centre[1] + blend_radius_m * sin(joint_az + blend_sweep * i / (blend_points - 1)),
        )
        for i in range(blend_points)
    ]
    beta = atan2(touch[0], touch[1] - offset_m)
    return [start] + arc + [joint] + blend, beta


def rotor_contour(
    geometry: RotaryGeometry,
    crank_angle_rad: float,
    mouth: MouthGeometry | None = None,
    od_points: int = 241,
    groove_points: int = 81,
) -> RotorContour:
    """Build the rotor boundary at one crank angle.

    Local coordinates are (transverse, axial) with the axial direction from
    the rotor center toward the groove center and the transverse direction
    positive toward the inlet.
    """

    mouth = MouthGeometry.default() if mouth is None else mouth
    mouth.validate(geometry)
    state = prescribed_state(geometry, crank_angle_rad)
    rotor_x, rotor_y = state.rotor_center_m
    groove_x, groove_y = state.cutout_center_m
    offset = geometry.cutout_offset_m
    axial = ((groove_x - rotor_x) / offset, (groove_y - rotor_y) / offset)
    transverse = (axial[1], -axial[0])

    def to_world(points: list[Point], mirror: bool = False) -> tuple[Point, ...]:
        sign = -1.0 if mirror else 1.0
        return tuple(
            (
                rotor_x + p[1] * axial[0] + sign * p[0] * transverse[0],
                rotor_y + p[1] * axial[1] + sign * p[0] * transverse[1],
            )
            for p in points
        )

    radius = geometry.rotor_radius_m
    flat_start = (
        radius * sin(mouth.flat_start_rad),
        radius * cos(mouth.flat_start_rad),
    )
    flat_end = (
        mouth.flat_end_radius_m * sin(mouth.flat_end_rad),
        mouth.flat_end_radius_m * cos(mouth.flat_end_rad),
    )
    run = (flat_end[0] - flat_start[0], flat_end[1] - flat_start[1])
    length = hypot(*run)
    inlet_lip, inlet_beta = _solve_lip(
        flat_end,
        (run[0] / length, run[1] / length),
        mouth.inlet_sweep_rad,
        offset,
        geometry.cutout_radius_m,
        mouth.lip_radius_m,
        mouth.blend_radius_m,
    )
    outlet_start = (
        radius * sin(mouth.outlet_tangent_rad),
        radius * cos(mouth.outlet_tangent_rad),
    )
    outlet_lip, outlet_beta = _solve_lip(
        outlet_start,
        (-cos(mouth.outlet_tangent_rad), sin(mouth.outlet_tangent_rad)),
        mouth.outlet_sweep_rad,
        offset,
        geometry.cutout_radius_m,
        mouth.lip_radius_m,
        mouth.blend_radius_m,
    )

    span = 2.0 * pi - mouth.flat_start_rad - mouth.outlet_tangent_rad
    od_local = [
        (
            radius * sin(mouth.flat_start_rad + span * i / (od_points - 1)),
            radius * cos(mouth.flat_start_rad + span * i / (od_points - 1)),
        )
        for i in range(od_points)
    ]
    groove_span = 2.0 * pi - inlet_beta - outlet_beta
    groove_local = [
        (
            geometry.cutout_radius_m * sin(inlet_beta + groove_span * i / (groove_points - 1)),
            offset
            + geometry.cutout_radius_m * cos(inlet_beta + groove_span * i / (groove_points - 1)),
        )
        for i in range(groove_points)
    ]

    od_arc = to_world(od_local)
    inlet_flat = to_world([flat_start, flat_end])
    mouth_path = (
        to_world(inlet_lip)
        + to_world(groove_local)
        + tuple(reversed(to_world(outlet_lip, mirror=True)))
    )
    material = od_arc + tuple(reversed(mouth_path))
    return RotorContour(
        od_arc=od_arc,
        inlet_flat=inlet_flat,
        mouth_path=mouth_path,
        material=material,
    )
