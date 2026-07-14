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

    Angle zero places the rotor at the top of its orbit. Positive angle advances
    clockwise. The circular-cutout center remains on the fixed vertical vane
    axis. An invisible rotor-fixed reference line prescribes the vane tip. The
    supplied
    tip distance is therefore a calibration at the top position, not a constant
    distance throughout the cycle.
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

    axis_x = horizontal_to_axis / geometry.cutout_offset_m
    axis_y = vertical_to_cutout / geometry.cutout_offset_m
    normal_x = -axis_y
    normal_y = axis_x
    reference_center_x = rotor_x + geometry.vane_tip_distance_at_top_m * axis_x
    reference_center_y = rotor_y + geometry.vane_tip_distance_at_top_m * axis_y
    if abs(normal_x) < 1.0e-12:
        raise ValueError("Vane reference line cannot intersect the vertical vane axis.")
    distance_along_reference = -reference_center_x / normal_x
    vane_tip_y = reference_center_y + distance_along_reference * normal_y
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


def distance(point_a: tuple[float, float], point_b: tuple[float, float]) -> float:
    """Return planar distance, primarily for validation and diagnostics."""

    return hypot(point_b[0] - point_a[0], point_b[1] - point_a[1])
