"""Port angular windows and the characteristic angles they define.

Ports were display markers until now: a single angle each, drawn on the bore
and used by nothing. This module turns them into angular windows and derives
the four characteristic angles the chamber pressure rule needs. See
PHYSICS.md section 3.4 for the supplied dimensions and their provenance.

All angles use the crank-angle convention of ``kinematics``: measured from
positive :math:`y` toward positive :math:`x`, which is clockwise on screen.
Because the geometry forces exact rotor-cylinder tangency, the contact point
sits at the bore angle equal to the crank angle, so a port window can be
compared against the crank angle without any conversion.

Do not mix this convention with ``PrescribedState.rotor_orientation_rad``,
which is a conventional ``atan2`` angle measured counterclockwise from
positive :math:`x`.
"""

from dataclasses import dataclass
from math import acos, isfinite, pi, radians, sqrt

from mochi.kinematics import RotaryGeometry

FULL_TURN_RAD = 2.0 * pi


@dataclass(frozen=True, slots=True)
class PortWindow:
    """A bore-wall opening spanning ``[start_rad, end_rad]`` clockwise.

    The window may wrap through zero, in which case ``end_rad`` is smaller
    than ``start_rad`` once both are normalized to ``[0, 2*pi)``.
    """

    start_rad: float
    end_rad: float

    @property
    def span_rad(self) -> float:
        """Clockwise angular width, always positive."""

        span = (self.end_rad - self.start_rad) % FULL_TURN_RAD
        return FULL_TURN_RAD if span == 0.0 else span

    @property
    def centre_rad(self) -> float:
        """Mid-angle of the window, normalized to ``[0, 2*pi)``."""

        return (self.start_rad + 0.5 * self.span_rad) % FULL_TURN_RAD


@dataclass(frozen=True, slots=True)
class CharacteristicAngles:
    """The four geometric characteristic angles of the cycle.

    Named after the supplied drawing symbols: ``suction_open_rad`` is
    :math:`\\varphi`, ``compression_start_rad`` is :math:`\\beta`, and the
    discharge port spans :math:`\\gamma` ending :math:`\\delta` before top
    dead center.
    """

    suction_open_rad: float
    compression_start_rad: float
    discharge_open_rad: float
    discharge_close_rad: float


def contains(window: PortWindow, angle_rad: float) -> bool:
    """Return True when ``angle_rad`` lies inside ``window``.

    Both bounds are inclusive and the wrap through zero is handled, so a
    window written as ``[350 deg, 10 deg]`` contains zero.
    """

    offset = (angle_rad - window.start_rad) % FULL_TURN_RAD
    return offset <= window.span_rad


def characteristic_angles(geometry: RotaryGeometry) -> CharacteristicAngles:
    """Derive the characteristic angles from the geometry's port timing."""

    geometry.validate()
    close_deg = 360.0 - geometry.recompression_angle_deg
    return CharacteristicAngles(
        suction_open_rad=radians(geometry.suction_seal_angle_deg),
        compression_start_rad=radians(geometry.compression_start_angle_deg),
        discharge_open_rad=radians(close_deg - geometry.discharge_port_span_deg),
        discharge_close_rad=radians(close_deg),
    )


def suction_window(geometry: RotaryGeometry) -> PortWindow:
    """Return the suction port as an angular window on the bore.

    The port opens at phi, its near edge, and its far edge is the
    compression-start angle beta where the chamber loses the suction port.
    phi is a machined port location, independent of the vane, so it is not
    expected to equal any vane-derived angle (PHYSICS.md section 3.4).
    """

    angles = characteristic_angles(geometry)
    return PortWindow(
        start_rad=angles.suction_open_rad,
        end_rad=angles.compression_start_rad,
    )


def discharge_window(geometry: RotaryGeometry) -> PortWindow:
    """Return the discharge port as an angular window on the bore."""

    angles = characteristic_angles(geometry)
    return PortWindow(
        start_rad=angles.discharge_open_rad,
        end_rad=angles.discharge_close_rad,
    )


def discharge_port_radius_m(geometry: RotaryGeometry) -> float:
    """Return the discharge port radius implied by its angular span.

    The port is treated as a circular opening whose diameter equals the arc
    its angular span subtends on the bore.
    """

    window = discharge_window(geometry)
    return 0.5 * geometry.cylinder_radius_m * window.span_rad


def port_open_area_m2(geometry: RotaryGeometry, crank_angle_rad: float) -> float:
    """Return the discharge port area still open to the compression chamber.

    The rotor-cylinder contact sweeps across the port and closes it from
    zero rather than all at once, so the exposed area is the circular
    segment of the port lying on the compression-chamber side of the
    contact. This is the geometric factor that a later orifice flow or
    valve model needs; no flow is computed here.
    """

    if not isfinite(crank_angle_rad):
        raise ValueError("Crank angle must be finite.")
    window = discharge_window(geometry)
    radius_m = discharge_port_radius_m(geometry)
    # The compression chamber runs clockwise from the contact to the vane at
    # 2*pi, so the port belongs to it while the contact angle is the smaller
    # of the two. Both angles are normalized first and compared directly:
    # wrapping the difference into (-pi, pi] would lose that ordering and
    # read the port as shut through most of compression.
    offset_rad = (crank_angle_rad % FULL_TURN_RAD) - window.centre_rad
    offset_m = geometry.cylinder_radius_m * offset_rad
    if offset_m >= radius_m:
        return 0.0
    if offset_m <= -radius_m:
        return pi * radius_m * radius_m
    return radius_m * radius_m * acos(offset_m / radius_m) - offset_m * sqrt(
        radius_m * radius_m - offset_m * offset_m
    )
