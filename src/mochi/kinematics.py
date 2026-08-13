"""Prescribed kinematics for the rotary-compressor test mechanism."""

from dataclasses import dataclass
from math import atan2, cos, hypot, isfinite, radians, sin, sqrt

MM = 1.0e-3


@dataclass(frozen=True, slots=True)
class RotaryGeometry:
    """Geometry and physical speed in SI units."""

    cylinder_id_m: float
    rotor_od_m: float
    eccentricity_m: float
    cutout_radius_m: float
    cutout_offset_m: float
    vane_width_m: float
    vane_tip_distance_at_top_m: float
    frequency_hz: float
    # Axial cylinder height supplied by the collaborators (2026-07-16); used
    # to turn chamber cross-section areas into volumes.
    cylinder_height_m: float = 21.0 * MM
    # Port timing supplied by the collaborators (2026-07-20), measured from
    # the vane centerline in the crank-angle sense (positive y toward
    # positive x). See PHYSICS.md section 3.4; these are independent supplied
    # numbers, not values derived from the vane width.
    suction_seal_angle_deg: float = 10.4
    compression_start_angle_deg: float = 27.7
    discharge_port_span_deg: float = 7.2
    recompression_angle_deg: float = 13.2
    # Stress-relief blend filleting each vane flank into the cylinder bore
    # (supplied 2026-07-23). See PHYSICS.md section 3.3.
    vane_cylinder_fillet_m: float = 2.1 * MM
    # Round on each corner of the fixed vane tip: the flank runs straight to
    # y_tip + this radius, then a quarter-circle to the tip flat. See
    # PHYSICS.md section 3.3.
    vane_tip_fillet_m: float = 1.5 * MM

    @classmethod
    def default(cls) -> "RotaryGeometry":
        """Return the dimensions supplied for the first test mechanism."""

        return cls(
            cylinder_id_m=77.0 * MM,
            rotor_od_m=68.0 * MM,
            eccentricity_m=4.5 * MM,
            cutout_radius_m=8.0 * MM,
            cutout_offset_m=25.0 * MM,
            vane_width_m=8.0 * MM,
            vane_tip_distance_at_top_m=9.0 * MM,
            frequency_hz=30.0,
            cylinder_height_m=21.0 * MM,
            suction_seal_angle_deg=10.4,
            compression_start_angle_deg=27.7,
            discharge_port_span_deg=7.2,
            recompression_angle_deg=13.2,
            vane_cylinder_fillet_m=2.1 * MM,
            vane_tip_fillet_m=1.5 * MM,
        )

    @property
    def cylinder_radius_m(self) -> float:
        return 0.5 * self.cylinder_id_m

    @property
    def rotor_radius_m(self) -> float:
        return 0.5 * self.rotor_od_m

    @property
    def radial_clearance_m(self) -> float:
        return self.cylinder_radius_m - self.rotor_radius_m

    @property
    def angular_speed_rad_s(self) -> float:
        return 2.0 * 3.141592653589793 * self.frequency_hz

    @property
    def speed_rpm(self) -> float:
        return 60.0 * self.frequency_hz

    def validate(self) -> None:
        """Raise ``ValueError`` when the prescribed geometry is impossible."""

        values = (
            self.cylinder_id_m,
            self.rotor_od_m,
            self.eccentricity_m,
            self.cutout_radius_m,
            self.cutout_offset_m,
            self.vane_width_m,
            self.vane_tip_distance_at_top_m,
            self.frequency_hz,
            self.cylinder_height_m,
            self.suction_seal_angle_deg,
            self.compression_start_angle_deg,
            self.discharge_port_span_deg,
            self.recompression_angle_deg,
            self.vane_cylinder_fillet_m,
            self.vane_tip_fillet_m,
        )
        if not all(isfinite(value) for value in values):
            raise ValueError("All geometry and speed values must be finite.")
        if self.cylinder_id_m <= 0.0 or self.rotor_od_m <= 0.0:
            raise ValueError("Cylinder ID and rotor OD must be positive.")
        if self.rotor_radius_m >= self.cylinder_radius_m:
            raise ValueError("Rotor OD must be smaller than cylinder ID.")
        if self.eccentricity_m < 0.0:
            raise ValueError("Eccentricity cannot be negative.")
        if self.eccentricity_m > self.radial_clearance_m + 1.0e-12:
            raise ValueError("Eccentricity makes the rotor intersect the cylinder.")
        if self.cutout_radius_m <= 0.0:
            raise ValueError("Cutout radius must be positive.")
        if self.cutout_offset_m <= self.eccentricity_m:
            raise ValueError("Cutout offset must exceed the eccentricity.")
        if self.cutout_offset_m + self.cutout_radius_m > self.rotor_radius_m:
            raise ValueError("Circular cutout must fit inside the rotor OD.")
        if self.vane_width_m <= 0.0:
            raise ValueError("Vane width must be positive.")
        if self.vane_width_m >= 2.0 * self.cutout_radius_m:
            raise ValueError("Vane must be narrower than the circular cutout.")
        if not 0.0 < self.vane_tip_distance_at_top_m < self.cutout_offset_m:
            raise ValueError(
                "Top-position vane-tip distance must lie between the rotor and cutout centers."
            )
        if self.frequency_hz <= 0.0:
            raise ValueError("Frequency must be positive.")
        if self.cylinder_height_m <= 0.0:
            raise ValueError("Cylinder height must be positive.")
        port_angles = (
            self.suction_seal_angle_deg,
            self.compression_start_angle_deg,
            self.discharge_port_span_deg,
            self.recompression_angle_deg,
        )
        if any(angle <= 0.0 for angle in port_angles):
            raise ValueError("Port timing angles must be positive.")
        if self.suction_seal_angle_deg >= self.compression_start_angle_deg:
            raise ValueError("Suction must open before compression starts.")
        discharge_open_deg = 360.0 - self.recompression_angle_deg - self.discharge_port_span_deg
        if discharge_open_deg <= self.compression_start_angle_deg:
            raise ValueError("Discharge port must open after compression starts.")
        if self.vane_cylinder_fillet_m <= 0.0:
            raise ValueError("Vane-cylinder blend radius must be positive.")
        if 2.0 * self.vane_cylinder_fillet_m >= self.cylinder_radius_m - 0.5 * self.vane_width_m:
            raise ValueError("Vane-cylinder blend radius is too large to meet the bore.")
        if self.vane_tip_fillet_m <= 0.0:
            raise ValueError("Vane-tip round radius must be positive.")
        if 2.0 * self.vane_tip_fillet_m >= self.vane_width_m:
            raise ValueError("Vane-tip round radius is too large to leave a tip flat.")


@dataclass(frozen=True, slots=True)
class PrescribedState:
    """Position and orientation for one crank angle."""

    crank_angle_rad: float
    rotor_center_m: tuple[float, float]
    cutout_center_m: tuple[float, float]
    rotor_orientation_rad: float
    vane_tip_m: tuple[float, float]
    vane_length_m: float


def prescribed_state(geometry: RotaryGeometry, crank_angle_rad: float) -> PrescribedState:
    """Evaluate the constrained rotor and circular-cutout position.

    Angle zero (top dead center) places the rotor at the top of its orbit.
    Positive angle advances clockwise. The circular-cutout center remains on
    the fixed vertical vane axis. The stepped vane is integral with the
    cylinder, so its tip is fixed in the global frame; the supplied tip
    distance is the rotor-center-to-tip distance at top dead center, which
    places the tip at y = eccentricity + tip distance at every crank angle.
    The apparent planar overlap with the rotor silhouette is resolved by the
    axial stepped structure (PHYSICS.md section 3.3).
    """

    geometry.validate()
    eccentricity = geometry.eccentricity_m
    rotor_x = eccentricity * sin(crank_angle_rad)
    rotor_y = eccentricity * cos(crank_angle_rad)

    horizontal_to_axis = -rotor_x
    vertical_to_cutout = sqrt(max(geometry.cutout_offset_m**2 - horizontal_to_axis**2, 0.0))
    cutout_x = 0.0
    cutout_y = rotor_y + vertical_to_cutout
    orientation = atan2(vertical_to_cutout, horizontal_to_axis)

    vane_tip_y = eccentricity + geometry.vane_tip_distance_at_top_m
    vane_tip = (0.0, vane_tip_y)
    vane_length = geometry.cylinder_radius_m - vane_tip_y

    return PrescribedState(
        crank_angle_rad=crank_angle_rad,
        rotor_center_m=(rotor_x, rotor_y),
        cutout_center_m=(cutout_x, cutout_y),
        rotor_orientation_rad=orientation,
        vane_tip_m=vane_tip,
        vane_length_m=vane_length,
    )


def port_position(
    cylinder_radius_m: float,
    angle_from_positive_y_deg: float,
) -> tuple[float, float]:
    """Return a cylinder-wall point for an angle measured from positive y."""

    angle = radians(angle_from_positive_y_deg)
    return (
        cylinder_radius_m * sin(angle),
        cylinder_radius_m * cos(angle),
    )


def vane_fillet_geometry(
    geometry: RotaryGeometry,
) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]]:
    """Right-side vane-root blend: ``(centre, flank_tangent, bore_tangent)``.

    The blend of radius ``vane_cylinder_fillet_m`` is tangent to the vane
    flank at ``x = +vane_width / 2`` and internally tangent to the bore. Its
    centre therefore sits one blend radius outside the flank and one blend
    radius inside the bore. The left-side blend is the mirror image across
    the vane centreline. Points are global metres.
    """

    fillet_m = geometry.vane_cylinder_fillet_m
    half_width_m = 0.5 * geometry.vane_width_m
    bore_m = geometry.cylinder_radius_m
    centre_x = half_width_m + fillet_m
    centre_y = sqrt((bore_m - fillet_m) ** 2 - centre_x**2)
    flank_tangent = (half_width_m, centre_y)
    scale = bore_m / (bore_m - fillet_m)
    bore_tangent = (centre_x * scale, centre_y * scale)
    return (centre_x, centre_y), flank_tangent, bore_tangent


def vane_tip_round(
    geometry: RotaryGeometry,
) -> tuple[float, float, tuple[float, float], tuple[float, float]]:
    """Right-corner vane-tip round: ``(tip_y, flank_top, centre, tip_tangent)``.

    Each corner of the fixed vane tip is rounded by ``vane_tip_fillet_m``: the
    flank at ``x = +vane_width / 2`` runs straight down to ``flank_top = y_tip +
    R``, then a quarter-circle of radius ``R`` curves in to the tip flat at
    ``y_tip``, tangent there at ``x = +(vane_width / 2 - R)``. Returns the tip
    ``y`` (from :func:`prescribed_state`, fixed in the global frame), the straight
    flank's lowest ``y`` (``flank_top``), the arc centre, and the tip-flat tangent
    point. The left corner is the mirror across ``x = 0``. Global metres.
    """

    radius_m = geometry.vane_tip_fillet_m
    half_width_m = 0.5 * geometry.vane_width_m
    tip_y = geometry.eccentricity_m + geometry.vane_tip_distance_at_top_m
    flank_top = tip_y + radius_m
    centre = (half_width_m - radius_m, flank_top)
    tip_tangent = (half_width_m - radius_m, tip_y)
    return tip_y, flank_top, centre, tip_tangent


def distance(point_a: tuple[float, float], point_b: tuple[float, float]) -> float:
    """Return planar distance, primarily for validation and diagnostics."""

    return hypot(point_b[0] - point_a[0], point_b[1] - point_a[1])
