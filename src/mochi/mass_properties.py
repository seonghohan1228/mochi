"""Solid-body mass and moment of inertia from the modelled cross-sections.

The rotor and swing-bush **material is not confirmed**, so the density is an
**editable assumption**: change ``ROTOR_DENSITY_KG_M3`` / ``BUSH_DENSITY_KG_M3`` (or
pass ``density_kg_m3``) and the mass and moment of inertia **rescale linearly** --
both are ``rho * integral dA * H`` and ``rho * integral r^2 dA * H``, so a density
edit changes them immediately, with no other input to touch.

The cross-sections are the *same solid regions the pressure/volume model already
classifies* (:mod:`mochi.chamber_volume`, :mod:`mochi.rotor_profile`): the rotor OD
disc minus the mouth cavity minus the central crank-pin bore, and the two swing-bush
pieces. They are integrated on the same raster as
:func:`mochi.chamber_volume.true_chamber_volumes` (mass/inertia are intrinsic body
properties, so they are evaluated once at ``theta = 0``).

This is the inertia data the D2 rotordynamics rung needs: the rotor/bush **lateral**
equations of motion and the orbital inertial load ``m_r e omega^2`` that the
quasi-static model dropped. See PHYSICS.md Section 4.8.
"""

from dataclasses import dataclass

import numpy as np

from mochi.bush_film import AXIAL_CLEARANCE_M
from mochi.chamber_volume import (
    DEFAULT_GRID_STEP_M,
    SwingBush,
    _mouth_cavity_polygon,
    _points_in_polygon,
)
from mochi.journal_bearing import JOURNAL_RADIUS_M
from mochi.kinematics import RotaryGeometry, prescribed_state

# EDITABLE assumed densities -- the rotor/bush material is not confirmed. Change
# these (or pass ``density_kg_m3``) and every mass/inertia rescales linearly.
ROTOR_DENSITY_KG_M3 = 7200.0  # grey cast iron (typical rolling-piston rotor)
BUSH_DENSITY_KG_M3 = 7850.0  # steel swing bush


@dataclass(frozen=True, slots=True)
class MassProperties:
    """Mass and polar moment of inertia of a solid body (SI units)."""

    mass_kg: float
    inertia_kg_m2: float  # polar (about +z) through the reference point named per function
    area_m2: float  # solid cross-section area
    height_m: float
    density_kg_m3: float
    grid_step_m: float


def _grid(geometry: RotaryGeometry, grid_step_m: float) -> tuple[np.ndarray, np.ndarray]:
    """Cell-centred square raster spanning the bore, matching `true_chamber_volumes`."""

    bore = geometry.cylinder_radius_m
    axis = np.arange(-bore, bore, grid_step_m) + 0.5 * grid_step_m
    return np.meshgrid(axis, axis)


def rotor_mass_properties(
    geometry: RotaryGeometry,
    *,
    density_kg_m3: float | None = None,
    bore_radius_m: float = JOURNAL_RADIUS_M,
    mouth=None,
    grid_step_m: float = DEFAULT_GRID_STEP_M,
) -> MassProperties:
    """Rotor solid mass and polar inertia **about the rotor centre O_r**.

    Solid = the OD disc (radius ``R_r`` about ``O_r``) minus the mouth cavity minus
    the central crank-pin bore (radius ``bore_radius_m``, default the journal radius
    ``r_j`` = 14.2 mm), at full cylinder height. ``density_kg_m3`` defaults to the
    editable :data:`ROTOR_DENSITY_KG_M3`; the polar inertia is the second area moment
    about ``O_r`` times ``rho * H`` (the axis the crank-pin journal and the D2 rotor
    tilt/whirl use).
    """

    geometry.validate()
    density = ROTOR_DENSITY_KG_M3 if density_kg_m3 is None else density_kg_m3
    if not (density > 0.0 and bore_radius_m > 0.0 and grid_step_m > 0.0):
        raise ValueError("Density, bore radius, and grid step must be positive.")
    if bore_radius_m >= geometry.rotor_radius_m:
        raise ValueError("Central bore must be smaller than the rotor outside diameter.")

    state = prescribed_state(geometry, 0.0)
    rotor_x, rotor_y = state.rotor_center_m
    grid_x, grid_y = _grid(geometry, grid_step_m)
    from_rotor = np.hypot(grid_x - rotor_x, grid_y - rotor_y)
    inside_od = from_rotor <= geometry.rotor_radius_m
    inside_bore = from_rotor <= bore_radius_m

    cavity = _mouth_cavity_polygon(geometry, 0.0, mouth)
    window = (
        (grid_x >= cavity[:, 0].min())
        & (grid_x <= cavity[:, 0].max())
        & (grid_y >= cavity[:, 1].min())
        & (grid_y <= cavity[:, 1].max())
    )
    in_cavity = np.zeros(grid_x.shape, dtype=bool)
    in_cavity[window] = _points_in_polygon(grid_x[window], grid_y[window], cavity)

    solid = inside_od & ~in_cavity & ~inside_bore
    cell_area = grid_step_m * grid_step_m
    height = geometry.cylinder_height_m
    area = float(solid.sum()) * cell_area
    r_squared = (grid_x - rotor_x) ** 2 + (grid_y - rotor_y) ** 2
    inertia = float((solid * r_squared).sum()) * cell_area * height * density
    return MassProperties(
        mass_kg=area * height * density,
        inertia_kg_m2=inertia,
        area_m2=area,
        height_m=height,
        density_kg_m3=density,
        grid_step_m=grid_step_m,
    )


def bush_mass_properties(
    geometry: RotaryGeometry,
    *,
    density_kg_m3: float | None = None,
    bush: SwingBush | None = None,
    grid_step_m: float = DEFAULT_GRID_STEP_M,
) -> MassProperties:
    """Swing-bush (both pieces) mass and polar inertia **about the groove centre**.

    Height = cylinder height minus twice the axial running clearance (the bush runs
    0.0085 mm short per side, PHYSICS.md 3.3). Both pieces via
    :meth:`mochi.chamber_volume.SwingBush.occupies`. ``density_kg_m3`` defaults to the
    editable :data:`BUSH_DENSITY_KG_M3`.
    """

    geometry.validate()
    density = BUSH_DENSITY_KG_M3 if density_kg_m3 is None else density_kg_m3
    if not (density > 0.0 and grid_step_m > 0.0):
        raise ValueError("Density and grid step must be positive.")
    bush = SwingBush() if bush is None else bush
    bush.validate(geometry)

    state = prescribed_state(geometry, 0.0)
    groove_x, groove_y = state.cutout_center_m
    grid_x, grid_y = _grid(geometry, grid_step_m)
    occupied = bush.occupies(grid_x, grid_y, state.cutout_center_m)

    cell_area = grid_step_m * grid_step_m
    height = geometry.cylinder_height_m - 2.0 * AXIAL_CLEARANCE_M
    area = float(occupied.sum()) * cell_area
    r_squared = (grid_x - groove_x) ** 2 + (grid_y - groove_y) ** 2
    inertia = float((occupied * r_squared).sum()) * cell_area * height * density
    return MassProperties(
        mass_kg=area * height * density,
        inertia_kg_m2=inertia,
        area_m2=area,
        height_m=height,
        density_kg_m3=density,
        grid_step_m=grid_step_m,
    )


@dataclass(frozen=True, slots=True)
class PieceMassProperties:
    """One swing-bush piece: mass, centre of mass, and inertia about the CM (SI)."""

    mass_kg: float
    inertia_about_cg_kg_m2: float
    cg_offset_m: tuple[float, float]  # piece centre of mass relative to the groove centre
    height_m: float
    density_kg_m3: float


def bush_piece_mass_properties(
    geometry: RotaryGeometry,
    side: float,
    *,
    density_kg_m3: float | None = None,
    bush: SwingBush | None = None,
    grid_step_m: float = DEFAULT_GRID_STEP_M,
) -> PieceMassProperties:
    """One bush piece's mass, CM, and polar inertia **about its own centre of mass**.

    ``side = +1`` selects the IN piece (``x > groove_x``), ``-1`` the OUT piece; the two
    pieces lie on opposite sides of the vane, so an ``x`` split isolates one. Returns
    the inertia about the piece CM (the axis the D2 piece-attitude EOM uses), the CM
    offset from the groove centre, and the mass -- the piece mass is editable /
    overridable by the confirmed 6.341 g (Section 4.11).
    """

    geometry.validate()
    if side not in (1.0, -1.0):
        raise ValueError("side must be +1 (IN) or -1 (OUT).")
    density = BUSH_DENSITY_KG_M3 if density_kg_m3 is None else density_kg_m3
    if not (density > 0.0 and grid_step_m > 0.0):
        raise ValueError("Density and grid step must be positive.")
    bush = SwingBush() if bush is None else bush
    bush.validate(geometry)

    state = prescribed_state(geometry, 0.0)
    groove_x, groove_y = state.cutout_center_m
    grid_x, grid_y = _grid(geometry, grid_step_m)
    occupied = bush.occupies(grid_x, grid_y, state.cutout_center_m) & (
        (grid_x > groove_x) if side > 0.0 else (grid_x < groove_x)
    )

    count = float(occupied.sum())
    if count <= 0.0:
        raise ValueError("No piece cells found; check the grid step and geometry.")
    cell_area = grid_step_m * grid_step_m
    height = geometry.cylinder_height_m - 2.0 * AXIAL_CLEARANCE_M
    cg_x = float((occupied * grid_x).sum()) / count
    cg_y = float((occupied * grid_y).sum()) / count
    r_squared = (grid_x - cg_x) ** 2 + (grid_y - cg_y) ** 2
    inertia = float((occupied * r_squared).sum()) * cell_area * height * density
    return PieceMassProperties(
        mass_kg=count * cell_area * height * density,
        inertia_about_cg_kg_m2=inertia,
        cg_offset_m=(cg_x - groove_x, cg_y - groove_y),
        height_m=height,
        density_kg_m3=density,
    )
