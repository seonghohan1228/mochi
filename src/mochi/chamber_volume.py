"""Chamber volumes from the true rotor contour, swing bush, and axial steps.

The circular-rotor approximation of :mod:`mochi.chambers` replaces the rotor
by its outside-diameter disc, which discards the gas space inside the rotor
mouth. That space is small against the swept volume but it is the whole of
the volume trapped when the discharge port closes, so the approximation
cannot support a recompression phase at all: it drives the trapped volume to
zero and the polytropic pressure to infinity.

This module integrates the real boundary instead: the rotor contour of
:mod:`mochi.rotor_profile`, the two swing-bush pieces, and the stepped vane.
See PHYSICS.md section 3.4.

Assumptions, all stated in PHYSICS.md section 3.4:

* The swing bush keeps a fixed attitude and only translates with the groove
  center. The rotor's relative tilt is not applied to the bush.
* The 0.010 mm bush films against the vane and the groove are treated as
  filled, so the bush seals the 4.0 MPa recess spaces away from the chamber
  gas and those spaces carry no chamber volume.
* The mouth runs through the full axial height; the vane is full thickness
  down to ``full_vane_depth_m`` from the bore and only its two ledges block
  the section below that.
"""

from dataclasses import dataclass
from math import atan2, cos, degrees, isfinite, pi, sin, sqrt

import numpy as np
from numpy.typing import NDArray

from mochi.kinematics import MM, RotaryGeometry, prescribed_state, vane_fillet_geometry
from mochi.ports import characteristic_angles
from mochi.rotor_profile import MouthGeometry, rotor_contour

FULL_TURN_DEG = 360.0
# Grid pitch of the area integration. The clearance volume changes by 0.3 %
# between this pitch and half of it, which sets the reported accuracy.
DEFAULT_GRID_STEP_M = 5.0e-5


@dataclass(frozen=True, slots=True)
class SwingBush:
    """Two half-moon pieces riding in the circular groove.

    Dimensions supplied 2026-07-20 (PHYSICS.md section 3.3). Each piece is
    the minor segment of a circle cut by a chord ``flat_offset_m`` from its
    own center, with both corners rounded by ``fillet_radius_m``. The pieces
    sit ``piece_shift_m`` toward their own side so the film against the vane
    and the film against the groove are both 0.010 mm.
    """

    piece_outer_radius_m: float = 7.970 * MM
    flat_offset_m: float = 3.990 * MM
    fillet_radius_m: float = 0.500 * MM
    piece_shift_m: float = 0.020 * MM

    def validate(self, geometry: RotaryGeometry) -> None:
        """Raise ``ValueError`` when the bush cannot sit in the groove."""

        values = (
            self.piece_outer_radius_m,
            self.flat_offset_m,
            self.fillet_radius_m,
            self.piece_shift_m,
        )
        if not all(isfinite(value) for value in values):
            raise ValueError("All swing-bush dimensions must be finite.")
        if self.piece_outer_radius_m <= 0.0 or self.fillet_radius_m <= 0.0:
            raise ValueError("Swing-bush radii must be positive.")
        if not 0.0 < self.flat_offset_m < self.piece_outer_radius_m:
            raise ValueError("Bush flat must lie inside its own piece radius.")
        if self.piece_outer_radius_m + self.piece_shift_m > geometry.cutout_radius_m:
            raise ValueError("Swing bush must fit inside the circular groove.")

    def half_arc_rad(self) -> float:
        """Return the half-angle of the cylindrical face, from the fillet."""

        return float(
            np.arccos(
                (self.flat_offset_m + self.fillet_radius_m)
                / (self.piece_outer_radius_m - self.fillet_radius_m)
            )
        )

    def occupies(
        self,
        x_m: NDArray[np.float64],
        y_m: NDArray[np.float64],
        groove_centre_m: tuple[float, float],
    ) -> NDArray[np.bool_]:
        """Return a mask of the points covered by either bush piece."""

        groove_x, groove_y = groove_centre_m
        covered = np.zeros(x_m.shape, dtype=bool)
        corner_x = self.flat_offset_m + self.fillet_radius_m
        corner_y = np.sqrt(
            (self.piece_outer_radius_m - self.fillet_radius_m) ** 2 - corner_x * corner_x
        )
        for side in (1.0, -1.0):
            # Each piece is shifted toward its own side of the vane.
            local_x = side * (x_m - (groove_x + side * self.piece_shift_m))
            local_y = y_m - groove_y
            piece = (np.hypot(local_x, local_y) <= self.piece_outer_radius_m) & (
                local_x >= self.flat_offset_m
            )
            fillet = np.minimum(
                np.hypot(local_x - corner_x, local_y - corner_y),
                np.hypot(local_x - corner_x, local_y + corner_y),
            )
            outside_fillet = (
                (local_x < corner_x)
                & (np.abs(local_y) > corner_y)
                & (fillet > self.fillet_radius_m)
            )
            covered |= piece & ~outside_fillet
        return covered


@dataclass(frozen=True, slots=True)
class AxialBands:
    """Axial structure of the stepped vane (PHYSICS.md section 3.3).

    The vane is full thickness from the bore down to ``full_vane_depth_m``
    and continues as two ledges of ``ledge_thickness_m`` each. Between the
    ledges the section is open to the chamber.
    """

    full_vane_depth_m: float = 15.4 * MM
    ledge_thickness_m: float = 2.4 * MM

    def validate(self, geometry: RotaryGeometry) -> None:
        """Raise ``ValueError`` when the axial steps do not fit the height."""

        if not isfinite(self.full_vane_depth_m) or not isfinite(self.ledge_thickness_m):
            raise ValueError("All axial band dimensions must be finite.")
        if self.full_vane_depth_m <= 0.0 or self.ledge_thickness_m <= 0.0:
            raise ValueError("Axial band dimensions must be positive.")
        if 2.0 * self.ledge_thickness_m >= geometry.cylinder_height_m:
            raise ValueError("Vane ledges must be thinner than the cylinder height.")

    def open_gap_m(self, geometry: RotaryGeometry) -> float:
        """Return the height left open between the two ledges."""

        return geometry.cylinder_height_m - 2.0 * self.ledge_thickness_m


@dataclass(frozen=True, slots=True)
class TrueChamberVolumes:
    """Chamber volumes at one crank angle from the real boundary."""

    crank_angle_rad: float
    suction_volume_m3: float
    discharge_volume_m3: float
    mouth_cavity_volume_m3: float
    grid_step_m: float


def _points_in_polygon(
    x_m: NDArray[np.float64],
    y_m: NDArray[np.float64],
    polygon: NDArray[np.float64],
) -> NDArray[np.bool_]:
    """Vectorized even-odd ray casting.

    NumPy is already a dependency; a polygon predicate does not justify
    adding a geometry library.
    """

    inside = np.zeros(x_m.shape, dtype=bool)
    start = polygon
    end = np.roll(polygon, -1, axis=0)
    for (x0, y0), (x1, y1) in zip(start, end, strict=True):
        if y0 == y1:
            continue
        straddles = (y0 > y_m) != (y1 > y_m)
        crossing_x = x0 + (y_m - y0) * (x1 - x0) / (y1 - y0)
        inside ^= straddles & (x_m < crossing_x)
    return inside


def _mouth_cavity_polygon(
    geometry: RotaryGeometry,
    crank_angle_rad: float,
    mouth: MouthGeometry | None,
) -> NDArray[np.float64]:
    """Return the gas cavity inside the rotor outside diameter.

    The rotor material is the outside-diameter disc minus this cavity, which
    is bounded by the short outside-diameter arc across the mouth, the inlet
    flat, and the mouth path through both lips and the groove.
    """

    contour = rotor_contour(geometry, crank_angle_rad, mouth=mouth)
    state = prescribed_state(geometry, crank_angle_rad)
    centre_x, centre_y = state.rotor_center_m
    azimuth = [atan2(y - centre_y, x - centre_x) for x, y in contour.od_arc]
    start_az = azimuth[-1]
    # The closing arc must keep turning the same way as the stored arc,
    # otherwise it retraces it and the polygon becomes the whole disc.
    step = (azimuth[1] - azimuth[0] + pi) % (2.0 * pi) - pi
    if step < 0.0:
        sweep = -((start_az - azimuth[0]) % (2.0 * pi))
    else:
        sweep = (azimuth[0] - start_az) % (2.0 * pi)
    radius = geometry.rotor_radius_m
    steps = 61
    short_arc = [
        (
            centre_x + radius * cos(start_az + sweep * index / (steps - 1)),
            centre_y + radius * sin(start_az + sweep * index / (steps - 1)),
        )
        for index in range(steps)
    ]
    cavity = short_arc + [contour.inlet_flat[1]] + list(contour.mouth_path)
    return np.asarray(cavity, dtype=float)


def true_chamber_volumes(
    geometry: RotaryGeometry,
    crank_angle_rad: float,
    bush: SwingBush | None = None,
    bands: AxialBands | None = None,
    mouth: MouthGeometry | None = None,
    grid_step_m: float = DEFAULT_GRID_STEP_M,
) -> TrueChamberVolumes:
    """Integrate the suction and discharge gas volumes at one crank angle.

    The split follows :mod:`mochi.chambers`: the suction chamber runs
    clockwise from the vane to the rotor-cylinder contact and the discharge
    chamber continues clockwise back to the vane. Unlike the circular-rotor
    approximation this stays finite through the whole revolution, because
    the mouth cavity never closes.
    """

    if not isfinite(crank_angle_rad):
        raise ValueError("Crank angle must be finite.")
    if not isfinite(grid_step_m) or grid_step_m <= 0.0:
        raise ValueError("Grid step must be positive and finite.")
    geometry.validate()
    bush = SwingBush() if bush is None else bush
    bands = AxialBands() if bands is None else bands
    bush.validate(geometry)
    bands.validate(geometry)

    state = prescribed_state(geometry, crank_angle_rad)
    rotor_x, rotor_y = state.rotor_center_m
    bore_radius = geometry.cylinder_radius_m
    axis = np.arange(-bore_radius, bore_radius, grid_step_m) + 0.5 * grid_step_m
    grid_x, grid_y = np.meshgrid(axis, axis)

    inside_bore = grid_x * grid_x + grid_y * grid_y <= bore_radius * bore_radius
    from_rotor = np.hypot(grid_x - rotor_x, grid_y - rotor_y)
    inside_rotor_od = from_rotor <= geometry.rotor_radius_m

    # The cavity polygon only spans the mouth, so the ray casting runs on a
    # small window instead of the whole bore.
    cavity = _mouth_cavity_polygon(geometry, crank_angle_rad, mouth)
    window = (
        (grid_x >= cavity[:, 0].min())
        & (grid_x <= cavity[:, 0].max())
        & (grid_y >= cavity[:, 1].min())
        & (grid_y <= cavity[:, 1].max())
    )
    in_cavity = np.zeros(grid_x.shape, dtype=bool)
    in_cavity[window] = _points_in_polygon(grid_x[window], grid_y[window], cavity)

    material = inside_rotor_od & ~in_cavity
    occupied_by_bush = bush.occupies(grid_x, grid_y, state.cutout_center_m)
    free = inside_bore & ~material & ~occupied_by_bush

    vane_column = (np.abs(grid_x) <= 0.5 * geometry.vane_width_m) & (grid_y >= state.vane_tip_m[1])
    full_vane = vane_column & (grid_y >= bore_radius - bands.full_vane_depth_m)
    height = np.where(
        vane_column & ~full_vane,
        bands.open_gap_m(geometry),
        geometry.cylinder_height_m,
    )
    height[full_vane] = 0.0
    # The vane-root blends are full-height solid on each side of the vane,
    # so they carry no gas: the pocket between the vane flank, the bore, and
    # the blend arc becomes material.
    (centre_x, centre_y), _, (bore_tangent_x, _) = vane_fillet_geometry(geometry)
    fillet_radius = geometry.vane_cylinder_fillet_m
    half_width = 0.5 * geometry.vane_width_m
    corner_y = sqrt(bore_radius**2 - half_width**2)
    for side in (1.0, -1.0):
        local_x = side * grid_x
        pocket = (
            (local_x >= half_width)
            & (local_x <= bore_tangent_x)
            & (grid_y >= centre_y)
            & (grid_y <= corner_y)
            & inside_bore
            & ((local_x - centre_x) ** 2 + (grid_y - centre_y) ** 2 >= fillet_radius**2)
        )
        height[pocket] = 0.0

    bore_angle_deg = np.degrees(np.arctan2(grid_x, grid_y)) % FULL_TURN_DEG
    contact_deg = degrees(crank_angle_rad) % FULL_TURN_DEG
    cell_area = grid_step_m * grid_step_m
    suction = free & (bore_angle_deg < contact_deg)
    discharge = free & ~suction
    return TrueChamberVolumes(
        crank_angle_rad=crank_angle_rad,
        suction_volume_m3=float((suction * height).sum() * cell_area),
        discharge_volume_m3=float((discharge * height).sum() * cell_area),
        mouth_cavity_volume_m3=float(((free & inside_rotor_od) * height).sum() * cell_area),
        grid_step_m=grid_step_m,
    )


def clearance_volume_m3(
    geometry: RotaryGeometry,
    bush: SwingBush | None = None,
    bands: AxialBands | None = None,
    mouth: MouthGeometry | None = None,
    grid_step_m: float = DEFAULT_GRID_STEP_M,
) -> float:
    """Return the volume trapped when the discharge port closes.

    This is the clearance volume the recompression phase starts from. It is
    a derived quantity, not a supplied one: the geometry fixes it.
    """

    angles = characteristic_angles(geometry)
    volumes = true_chamber_volumes(
        geometry,
        angles.discharge_close_rad,
        bush=bush,
        bands=bands,
        mouth=mouth,
        grid_step_m=grid_step_m,
    )
    return volumes.discharge_volume_m3


def mouth_cavity_area_m2(
    geometry: RotaryGeometry,
    crank_angle_rad: float,
    bush: SwingBush | None = None,
    mouth: MouthGeometry | None = None,
    grid_step_m: float = DEFAULT_GRID_STEP_M,
) -> float:
    """Return the free cross-section area inside the rotor mouth.

    PHYSICS.md section 3.3 recorded this as roughly 30 mm^2 while it was a
    neglected volume; this function is the independent measurement of it.
    """

    volumes = true_chamber_volumes(
        geometry,
        crank_angle_rad,
        bush=bush,
        mouth=mouth,
        grid_step_m=grid_step_m,
    )
    return volumes.mouth_cavity_volume_m3 / geometry.cylinder_height_m
