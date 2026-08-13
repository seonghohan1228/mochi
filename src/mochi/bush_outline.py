"""Shared 2-D outlines (in millimetres) for the swing-bush pieces and the vane tongue.

Both the static ``results/`` figures (:mod:`scripts.generate_results`) and the coupled-orbit
animation (:mod:`mochi.bush_gui`) need the *same* bush-piece and vane profiles as plain
``(xs, ys)`` millimetre point lists. Keeping the construction in one place stops the static
figures and the animation from silently diverging when the geometry is tweaked. These mirror
the canvas-space constructions drawn by :mod:`mochi.gui` (``_draw_swing_bush`` /
``_draw_cylinder_vane_outline``), which render straight onto a Tk canvas rather than returning
point lists.
"""

from __future__ import annotations

from math import atan2, cos, pi, sin, sqrt

from mochi.chamber_volume import SwingBush
from mochi.kinematics import MM, RotaryGeometry, vane_fillet_geometry


def bush_piece_outline_mm(
    centre_x_m: float, groove_y_m: float, side: float
) -> tuple[list[float], list[float]]:
    """One swing-bush piece outline in mm (mirrors ``gui._draw_swing_bush``)."""

    bush = SwingBush()
    radius = bush.piece_outer_radius_m
    fillet = bush.fillet_radius_m
    half_arc = bush.half_arc_rad()
    corner_x = bush.flat_offset_m + fillet
    corner_y = sqrt((radius - fillet) ** 2 - corner_x**2)
    fillet_start = atan2(radius * sin(half_arc) - corner_y, radius * cos(half_arc) - corner_x)
    xs: list[float] = []
    ys: list[float] = []

    def push(x_m: float, y_m: float) -> None:
        xs.append(x_m / MM)
        ys.append(y_m / MM)

    arc_steps, fillet_steps = 41, 8
    for index in range(arc_steps):
        angle = -half_arc + 2.0 * half_arc * index / (arc_steps - 1)
        push(centre_x_m + side * radius * cos(angle), groove_y_m + radius * sin(angle))
    for index in range(1, fillet_steps + 1):
        angle = fillet_start + (pi - fillet_start) * index / fillet_steps
        push(
            centre_x_m + side * (corner_x + fillet * cos(angle)),
            groove_y_m + corner_y + fillet * sin(angle),
        )
    push(centre_x_m + side * bush.flat_offset_m, groove_y_m + corner_y)
    push(centre_x_m + side * bush.flat_offset_m, groove_y_m - corner_y)
    for index in range(1, fillet_steps + 1):
        angle = pi + (pi - fillet_start) * index / fillet_steps
        push(
            centre_x_m + side * (corner_x + fillet * cos(angle)),
            groove_y_m - corner_y + fillet * sin(angle),
        )
    return xs, ys


def vane_outline_mm(
    geometry: RotaryGeometry, vane_tip_y_m: float, *, close_bore_arc: bool = True
) -> tuple[list[float], list[float]]:
    """Vane tongue outline in mm (mirrors ``gui._draw_cylinder_vane_outline``).

    ``close_bore_arc`` traces the short cylinder-bore arc across the top of the two root
    fillets so the filled polygon follows the bore (the static-figure convention). The
    zoomed-out animation leaves it off so the tongue closes with a straight chord, matching
    the profile that viewer has always drawn.
    """

    half_width = 0.5 * geometry.vane_width_m
    fillet = geometry.vane_cylinder_fillet_m
    (centre_x, centre_y), _, (bore_x, bore_y) = vane_fillet_geometry(geometry)
    blend_points = 16
    xs: list[float] = []
    ys: list[float] = []

    def push(x_m: float, y_m: float) -> None:
        xs.append(x_m / MM)
        ys.append(y_m / MM)

    # Right blend, bore tangent down to the flank tangent.
    right_start = atan2(bore_y - centre_y, bore_x - centre_x)
    for index in range(blend_points + 1):
        angle = right_start + (pi - right_start) * index / blend_points
        push(centre_x + fillet * cos(angle), centre_y + fillet * sin(angle))
    # Down the right flank to the tip round, across the shortened tip flat, up the left flank.
    tip_radius = geometry.vane_tip_fillet_m
    flank_top = vane_tip_y_m + tip_radius
    inner = half_width - tip_radius
    arc_steps = 8
    push(half_width, flank_top)
    for index in range(1, arc_steps + 1):
        angle = -0.5 * pi * index / arc_steps
        push(inner + tip_radius * cos(angle), flank_top + tip_radius * sin(angle))
    for index in range(arc_steps + 1):
        angle = -0.5 * pi - 0.5 * pi * index / arc_steps
        push(-inner + tip_radius * cos(angle), flank_top + tip_radius * sin(angle))
    # Left blend, flank tangent up to the bore tangent (index 0 lands on (-half_width, centre_y),
    # the flank tangent point). Same traversal order as gui._draw_cylinder_vane_outline.
    left_end = atan2(bore_y - centre_y, centre_x - bore_x)
    for index in range(blend_points + 1):
        angle = left_end * index / blend_points
        push(-centre_x + fillet * cos(angle), centre_y + fillet * sin(angle))
    if close_bore_arc:
        # Short bore arc across the top, left tangent back to the right tangent.
        bore = geometry.cylinder_radius_m
        right_bore_az = atan2(bore_y, bore_x)
        left_bore_az = pi - right_bore_az
        top_steps = 24
        for index in range(top_steps + 1):
            angle = left_bore_az + (right_bore_az - left_bore_az) * index / top_steps
            push(bore * cos(angle), bore * sin(angle))
    return xs, ys
