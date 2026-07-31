"""Tkinter test GUI for prescribed rotary-compressor motion."""

from __future__ import annotations

import math
import time
import tkinter as tk
from collections.abc import Callable
from tkinter import font as tkfont
from tkinter import ttk

from mochi.chamber_volume import SwingBush
from mochi.chambers import (
    ChamberSeparationError,
    chamber_areas,
    chamber_pressures,
    in_seal_over_window,
    seal_over_half_angle_rad,
)
from mochi.kinematics import (
    MM,
    PrescribedState,
    RotaryGeometry,
    port_position,
    prescribed_state,
    vane_fillet_geometry,
)
from mochi.ports import (
    CharacteristicAngles,
    characteristic_angles,
    discharge_window,
    port_open_area_m2,
    suction_window,
)
from mochi.rotor_profile import rotor_contour

CanvasTransform = Callable[[float, float], tuple[float, float]]


# Character count of the chamber status line. Every branch is padded to it so
# the fixed-width status label keeps a constant requested width.
STATUS_WIDTH = 127
# Same reason for the port status line, which also changes length by branch.
PORT_STATUS_WIDTH = 118
FRAME_DELAY_MS = 16
DEFAULT_SLOW_FACTOR = 0.01
FULL_TURN_RAD = 2.0 * math.pi


def stop_angle_crossed(
    previous_angle_rad: float,
    advanced_angle_rad: float,
    target_angle_rad: float,
) -> bool:
    """Return True when one clockwise animation step passes the target angle.

    The target is compared modulo one revolution, so any frame step size and
    wrap-around at 360 degrees are handled. Starting exactly on the target does
    not trigger, so resuming from a stop completes a full revolution before
    stopping there again.
    """

    step = advanced_angle_rad - previous_angle_rad
    if step <= 0.0:
        return False
    remaining = (target_angle_rad - previous_angle_rad) % FULL_TURN_RAD
    if remaining == 0.0:
        remaining = FULL_TURN_RAD
    return remaining <= step


class RotaryCompressorApp:
    """Animate the supplied rotary-compressor geometry on a Tk canvas."""

    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("mochi - rotary compressor prescribed-motion test")
        self.root.minsize(980, 700)
        # An explicit size stops Tk from resizing the window whenever the
        # status text changes length, which would rescale the drawing on
        # every frame. The user can still resize the window freely.
        self.root.geometry("1180x780")

        self.geometry = RotaryGeometry.default()
        self.slow_factor = DEFAULT_SLOW_FACTOR
        self.crank_angle_rad = 0.0
        self.stop_angle_rad: float | None = None
        self.running = False
        self.last_tick_s = time.perf_counter()

        self.cylinder_id_var = tk.StringVar(value="77.0")
        self.rotor_od_var = tk.StringVar(value="68.0")
        self.eccentricity_var = tk.StringVar(value="4.5")
        self.frequency_var = tk.StringVar(value="30.0")
        self.slow_factor_var = tk.StringVar(value=f"{DEFAULT_SLOW_FACTOR:g}")
        self.cutout_radius_var = tk.StringVar(value="8.0")
        self.cutout_offset_var = tk.StringVar(value="25.0")
        self.vane_width_var = tk.StringVar(value="8.0")
        self.vane_tip_distance_at_top_var = tk.StringVar(value="9.0")
        self.suction_seal_var = tk.StringVar(value="10.4")
        self.compression_start_var = tk.StringVar(value="27.7")
        self.discharge_span_var = tk.StringVar(value="7.2")
        self.recompression_var = tk.StringVar(value="13.2")
        self.stop_angle_var = tk.StringVar(value="0.0")
        self.stop_enabled_var = tk.BooleanVar(value=False)
        self.lock_eccentricity_var = tk.BooleanVar(value=True)
        self.status_var = tk.StringVar()
        self.error_var = tk.StringVar()
        self.start_text_var = tk.StringVar(value="Start")

        self._build_layout()
        self.canvas.bind("<Configure>", lambda _event: self._redraw())
        self.root.protocol("WM_DELETE_WINDOW", self.root.destroy)
        self._redraw()
        self.root.after(FRAME_DELAY_MS, self._tick)

    def _build_layout(self) -> None:
        outer = ttk.Frame(self.root, padding=12)
        outer.pack(fill=tk.BOTH, expand=True)

        controls = ttk.LabelFrame(outer, text="Prescribed inputs", padding=12)
        controls.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 12))

        fields = (
            ("Cylinder ID (mm)", self.cylinder_id_var),
            ("Rotor OD (mm)", self.rotor_od_var),
            ("Eccentricity (mm)", self.eccentricity_var),
            ("Physical speed (Hz)", self.frequency_var),
            ("Display slow factor", self.slow_factor_var),
            ("Circular cutout radius (mm)", self.cutout_radius_var),
            ("Cutout center distance (mm)", self.cutout_offset_var),
            ("Vane width (mm)", self.vane_width_var),
            ("Vane tip distance at top (mm)", self.vane_tip_distance_at_top_var),
            ("Suction seal phi (deg)", self.suction_seal_var),
            ("Compression start beta (deg)", self.compression_start_var),
            ("Discharge port span gamma (deg)", self.discharge_span_var),
            ("Recompression delta (deg)", self.recompression_var),
            ("Stop at crank angle (deg)", self.stop_angle_var),
        )
        for row, (label, variable) in enumerate(fields):
            ttk.Label(controls, text=label).grid(row=row, column=0, sticky=tk.W, pady=3)
            ttk.Entry(controls, textvariable=variable, width=13).grid(
                row=row,
                column=1,
                sticky=tk.EW,
                padx=(8, 0),
                pady=3,
            )

        lock = ttk.Checkbutton(
            controls,
            text="Lock e = (cylinder ID - rotor OD) / 2",
            variable=self.lock_eccentricity_var,
            command=self._sync_eccentricity,
        )
        lock.grid(row=len(fields), column=0, columnspan=2, sticky=tk.W, pady=(8, 4))

        stop = ttk.Checkbutton(
            controls,
            text="Stop at the crank angle above",
            variable=self.stop_enabled_var,
        )
        stop.grid(row=len(fields) + 1, column=0, columnspan=2, sticky=tk.W, pady=(0, 4))

        button_row = ttk.Frame(controls)
        button_row.grid(row=len(fields) + 2, column=0, columnspan=2, sticky=tk.EW, pady=(8, 4))
        ttk.Button(button_row, text="Apply", command=self._apply_inputs).pack(
            side=tk.LEFT,
            expand=True,
            fill=tk.X,
        )
        ttk.Button(
            button_row,
            textvariable=self.start_text_var,
            command=self._toggle_running,
        ).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=5)
        ttk.Button(button_row, text="Reset", command=self._reset).pack(
            side=tk.LEFT,
            expand=True,
            fill=tk.X,
        )

        ttk.Separator(controls).grid(
            row=len(fields) + 3,
            column=0,
            columnspan=2,
            sticky=tk.EW,
            pady=10,
        )
        ttk.Label(
            controls,
            text=(
                "Ports are angular windows on the\n"
                "bore, measured from +y toward +x:\n"
                "suction opens at phi and closes\n"
                "at beta, discharge spans gamma\n"
                "and shuts delta before top.\n\n"
                "Slow factor 0.01 displays\n"
                "motion at 1/100 real speed.\n\n"
                "Keep 'Lock e' checked: with\n"
                "e = (ID - OD) / 2 the rotor seals\n"
                "on the cylinder and the chamber\n"
                "split stays defined.\n\n"
                "The chamber line uses the\n"
                "circular-rotor approximation.\n"
                "True-geometry volumes including\n"
                "the mouth and bush are too slow\n"
                "to draw per frame. Pressures are\n"
                "absolute."
            ),
            justify=tk.LEFT,
        ).grid(row=len(fields) + 4, column=0, columnspan=2, sticky=tk.W)

        ttk.Label(
            controls,
            textvariable=self.error_var,
            foreground="#b42318",
            wraplength=265,
            justify=tk.LEFT,
        ).grid(row=len(fields) + 5, column=0, columnspan=2, sticky=tk.W, pady=(12, 0))

        view = ttk.Frame(outer)
        view.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.canvas = tk.Canvas(
            view,
            background="white",
            highlightthickness=1,
            highlightbackground="#c7ccd4",
        )
        self.canvas.pack(fill=tk.BOTH, expand=True)
        status_font = tkfont.nametofont("TkFixedFont").copy()
        status_font.configure(size=8)
        ttk.Label(
            view,
            textvariable=self.status_var,
            anchor=tk.W,
            justify=tk.LEFT,
            padding=(4, 8),
            font=status_font,
        ).pack(fill=tk.X)

    def _sync_eccentricity(self) -> None:
        if not self.lock_eccentricity_var.get():
            return
        try:
            cylinder_id = float(self.cylinder_id_var.get())
            rotor_od = float(self.rotor_od_var.get())
        except ValueError:
            return
        self.eccentricity_var.set(f"{0.5 * (cylinder_id - rotor_od):.6g}")

    def _read_inputs(self) -> tuple[RotaryGeometry, float, float | None]:
        self._sync_eccentricity()
        try:
            geometry = RotaryGeometry(
                cylinder_id_m=float(self.cylinder_id_var.get()) * MM,
                rotor_od_m=float(self.rotor_od_var.get()) * MM,
                eccentricity_m=float(self.eccentricity_var.get()) * MM,
                frequency_hz=float(self.frequency_var.get()),
                cutout_radius_m=float(self.cutout_radius_var.get()) * MM,
                cutout_offset_m=float(self.cutout_offset_var.get()) * MM,
                vane_width_m=float(self.vane_width_var.get()) * MM,
                vane_tip_distance_at_top_m=float(self.vane_tip_distance_at_top_var.get()) * MM,
                suction_seal_angle_deg=float(self.suction_seal_var.get()),
                compression_start_angle_deg=float(self.compression_start_var.get()),
                discharge_port_span_deg=float(self.discharge_span_var.get()),
                recompression_angle_deg=float(self.recompression_var.get()),
            )
            slow_factor = float(self.slow_factor_var.get())
            stop_angle_deg = float(self.stop_angle_var.get())
        except ValueError as error:
            raise ValueError("All inputs must be valid numbers.") from error
        geometry.validate()
        if not math.isfinite(slow_factor) or not 0.0 < slow_factor <= 1.0:
            raise ValueError("Display slow factor must be greater than 0 and at most 1.")
        stop_angle_rad: float | None = None
        if self.stop_enabled_var.get():
            if not math.isfinite(stop_angle_deg):
                raise ValueError("Stop crank angle must be a finite number of degrees.")
            stop_angle_rad = math.radians(stop_angle_deg) % FULL_TURN_RAD
        return geometry, slow_factor, stop_angle_rad

    def _apply_inputs(self) -> None:
        try:
            geometry, slow_factor, stop_angle_rad = self._read_inputs()
        except ValueError as error:
            self.error_var.set(str(error))
            return
        self.geometry = geometry
        self.slow_factor = slow_factor
        self.stop_angle_rad = stop_angle_rad
        self.error_var.set("")
        self._redraw()

    def _toggle_running(self) -> None:
        if not self.running:
            self._apply_inputs()
            if self.error_var.get():
                return
        self.running = not self.running
        self.start_text_var.set("Pause" if self.running else "Start")
        self.last_tick_s = time.perf_counter()

    def _reset(self) -> None:
        self.running = False
        self.start_text_var.set("Start")
        self.crank_angle_rad = 0.0
        self.last_tick_s = time.perf_counter()
        self._redraw()

    def _tick(self) -> None:
        now_s = time.perf_counter()
        elapsed_s = now_s - self.last_tick_s
        self.last_tick_s = now_s
        if self.running:
            previous = self.crank_angle_rad
            advanced = previous + self.geometry.angular_speed_rad_s * self.slow_factor * elapsed_s
            if self.stop_angle_rad is not None and stop_angle_crossed(
                previous, advanced, self.stop_angle_rad
            ):
                advanced = self.stop_angle_rad
                self.running = False
                self.start_text_var.set("Start")
            self.crank_angle_rad = advanced % FULL_TURN_RAD
            self._redraw()
        self.root.after(FRAME_DELAY_MS, self._tick)

    def _redraw(self) -> None:
        if not hasattr(self, "canvas"):
            return
        self.canvas.delete("all")
        width = max(self.canvas.winfo_width(), 640)
        height = max(self.canvas.winfo_height(), 560)
        geometry = self.geometry
        state = prescribed_state(geometry, self.crank_angle_rad)

        cylinder_radius = geometry.cylinder_radius_m
        top = cylinder_radius + 6.0 * MM
        bottom = -cylinder_radius - 6.0 * MM
        half_width = cylinder_radius + 6.0 * MM
        margin = 35.0
        scale = min(
            (width - 2.0 * margin) / (2.0 * half_width),
            (height - 2.0 * margin) / (top - bottom),
        )
        origin_x = 0.5 * width
        origin_y = margin + top * scale

        def point(x_m: float, y_m: float) -> tuple[float, float]:
            return origin_x + x_m * scale, origin_y - y_m * scale

        self._draw_rotor_outline(point, geometry, state)
        self._draw_swing_bush(point, state)
        self._draw_cylinder_vane_outline(point, geometry, state.vane_tip_m)
        self._draw_center_guides(point, geometry, state)

        # Ports are angular windows on the bore, drawn as arcs outside the
        # cylinder outline so they annotate it without masking it.
        angles = characteristic_angles(geometry)
        suction = suction_window(geometry)
        discharge = discharge_window(geometry)
        self._draw_port_arc(
            point,
            geometry,
            math.degrees(suction.start_rad),
            math.degrees(suction.end_rad),
            "IN",
        )
        self._draw_port_arc(
            point,
            geometry,
            math.degrees(discharge.start_rad),
            math.degrees(discharge.end_rad),
            "OUT",
        )

        self._set_status(state, angles)

    def _draw_center_guides(
        self,
        point: CanvasTransform,
        geometry: RotaryGeometry,
        state: PrescribedState,
    ) -> None:
        """Draw the cylinder center and dashed rotor-center orbit."""

        eccentricity = geometry.eccentricity_m
        locus_top_left = point(-eccentricity, eccentricity)
        locus_bottom_right = point(eccentricity, -eccentricity)
        guide_color = "#5c6068"
        self.canvas.create_oval(
            *locus_top_left,
            *locus_bottom_right,
            fill="",
            outline=guide_color,
            width=1,
            dash=(5, 4),
        )

        for center, label in (((0.0, 0.0), "C"), (state.rotor_center_m, "R")):
            center_x, center_y = point(*center)
            marker_radius_px = 4
            self.canvas.create_line(
                center_x - marker_radius_px,
                center_y,
                center_x + marker_radius_px,
                center_y,
                fill=guide_color,
                width=2,
            )
            self.canvas.create_line(
                center_x,
                center_y - marker_radius_px,
                center_x,
                center_y + marker_radius_px,
                fill=guide_color,
                width=2,
            )
            self.canvas.create_text(
                center_x + 6,
                center_y - 6,
                text=label,
                anchor=tk.SW,
                fill=guide_color,
                font=("TkDefaultFont", 9, "bold"),
            )

    def _draw_port_arc(
        self,
        point: CanvasTransform,
        geometry: RotaryGeometry,
        start_deg: float,
        end_deg: float,
        label: str,
    ) -> None:
        """Stroke one port as an arc on the bore with end ticks and a label."""

        radius = geometry.cylinder_radius_m + 2.0 * MM
        span_deg = (end_deg - start_deg) % 360.0
        steps = max(int(span_deg), 2) + 1
        arc: list[float] = []
        for index in range(steps):
            angle_deg = start_deg + span_deg * index / (steps - 1)
            arc.extend(point(*port_position(radius, angle_deg)))
        self.canvas.create_line(arc, fill="black", width=3, capstyle=tk.ROUND)
        for angle_deg in (start_deg, end_deg):
            inner = port_position(geometry.cylinder_radius_m - 1.0 * MM, angle_deg)
            outer = port_position(radius + 1.5 * MM, angle_deg)
            self.canvas.create_line(*point(*inner), *point(*outer), fill="black", width=2)
        label_point = port_position(
            geometry.cylinder_radius_m + 7.0 * MM, start_deg + 0.5 * span_deg
        )
        self.canvas.create_text(
            *point(*label_point), text=label, fill="black", font=("TkDefaultFont", 10, "bold")
        )

    def _draw_swing_bush(self, point: CanvasTransform, state: PrescribedState) -> None:
        """Draw both swing-bush pieces at the groove center.

        Attitude is fixed and only the position follows the groove, which is
        the assumption the volume model uses (PHYSICS.md section 3.4).
        """

        bush = SwingBush()
        groove_x, groove_y = state.cutout_center_m
        radius = bush.piece_outer_radius_m
        fillet = bush.fillet_radius_m
        half_arc = bush.half_arc_rad()
        # R0.5 fillet corner, tangent to the flat on one side and the outer
        # face on the other -- the same construction as SwingBush.occupies.
        corner_x = bush.flat_offset_m + fillet
        corner_y = math.sqrt((radius - fillet) ** 2 - corner_x**2)
        fillet_start = math.atan2(
            radius * math.sin(half_arc) - corner_y,
            radius * math.cos(half_arc) - corner_x,
        )
        arc_steps = 41
        fillet_steps = 8
        for side in (1.0, -1.0):
            centre_x = groove_x + side * bush.piece_shift_m
            outline: list[float] = []
            # Outer cylindrical face, from the bottom fillet up to the top one.
            for index in range(arc_steps):
                angle = -half_arc + 2.0 * half_arc * index / (arc_steps - 1)
                outline.extend(
                    point(
                        centre_x + side * radius * math.cos(angle),
                        groove_y + radius * math.sin(angle),
                    )
                )
            # Top fillet, rounding the outer face into the flat.
            for index in range(1, fillet_steps + 1):
                angle = fillet_start + (math.pi - fillet_start) * index / fillet_steps
                outline.extend(
                    point(
                        centre_x + side * (corner_x + fillet * math.cos(angle)),
                        groove_y + corner_y + fillet * math.sin(angle),
                    )
                )
            # Flat face, top tangent point to bottom tangent point.
            outline.extend(point(centre_x + side * bush.flat_offset_m, groove_y + corner_y))
            outline.extend(point(centre_x + side * bush.flat_offset_m, groove_y - corner_y))
            # Bottom fillet, rounding the flat back into the outer face.
            for index in range(1, fillet_steps + 1):
                angle = math.pi + (math.pi - fillet_start) * index / fillet_steps
                outline.extend(
                    point(
                        centre_x + side * (corner_x + fillet * math.cos(angle)),
                        groove_y - corner_y + fillet * math.sin(angle),
                    )
                )
            self.canvas.create_polygon(
                outline,
                fill="#b9c0cb",
                outline="black",
                width=1,
            )

    def _draw_rotor_outline(
        self,
        point: CanvasTransform,
        geometry: RotaryGeometry,
        state: PrescribedState,
    ) -> None:
        """Draw the confirmed rotor contour (PHYSICS.md section 3.3).

        The mouth is asymmetric: the inlet side starts from a straight cut in
        the outside diameter, the outlet side from an OD tangent. Only the real
        material edges are stroked, so the mouth stays open to the chamber.
        """

        contour = rotor_contour(geometry, state.crank_angle_rad)

        def canvas_points(points: tuple[tuple[float, float], ...]) -> list[float]:
            flat: list[float] = []
            for world_point in points:
                flat.extend(point(*world_point))
            return flat

        self.canvas.create_polygon(
            canvas_points(contour.material),
            fill="#d9d9d9",
            outline="",
        )
        for edge in (contour.od_arc, contour.inlet_flat, contour.mouth_path):
            self.canvas.create_line(
                canvas_points(edge),
                fill="black",
                width=3,
                joinstyle=tk.ROUND,
                capstyle=tk.ROUND,
            )

    def _draw_cylinder_vane_outline(
        self,
        point: CanvasTransform,
        geometry: RotaryGeometry,
        vane_tip: tuple[float, float],
    ) -> None:
        """Draw the cylinder and changing-length vane as one transparent outline.

        Each vane flank meets the bore through the confirmed blend radius
        (PHYSICS.md section 3.3), so the two roots are rounded rather than
        sharp.
        """

        radius = geometry.cylinder_radius_m
        half_vane_width = 0.5 * geometry.vane_width_m
        fillet = geometry.vane_cylinder_fillet_m
        (centre_x, centre_y), _, (bore_x, bore_y) = vane_fillet_geometry(geometry)
        right_angle = math.atan2(bore_y, bore_x)
        left_angle = math.pi - right_angle
        blend_points = 16

        outline: list[tuple[float, float]] = []
        # Bore arc, from the left blend's bore tangent the long way round to
        # the right blend's bore tangent.
        for index in range(181):
            angle = left_angle + (2.0 * math.pi + right_angle - left_angle) * index / 180
            outline.append((radius * math.cos(angle), radius * math.sin(angle)))
        # Right blend, bore tangent to flank tangent.
        right_start = math.atan2(bore_y - centre_y, bore_x - centre_x)
        for index in range(1, blend_points + 1):
            angle = right_start + (math.pi - right_start) * index / blend_points
            outline.append(
                (centre_x + fillet * math.cos(angle), centre_y + fillet * math.sin(angle))
            )
        # Down the right flank to the R1.5 tip round, across the shortened tip
        # flat, up the left flank (round leaves a tip flat of width w - 2R).
        tip_radius = geometry.vane_tip_fillet_m
        flank_top = vane_tip[1] + tip_radius
        inner = half_vane_width - tip_radius
        arc_steps = 8
        outline.append((half_vane_width, flank_top))
        for index in range(1, arc_steps + 1):
            angle = -0.5 * math.pi * index / arc_steps
            outline.append(
                (inner + tip_radius * math.cos(angle), flank_top + tip_radius * math.sin(angle))
            )
        for index in range(arc_steps + 1):
            angle = -0.5 * math.pi - 0.5 * math.pi * index / arc_steps
            outline.append(
                (-inner + tip_radius * math.cos(angle), flank_top + tip_radius * math.sin(angle))
            )
        outline.append((-half_vane_width, centre_y))
        # Left blend, flank tangent to bore tangent, mirrored.
        left_end = math.atan2(bore_y - centre_y, centre_x - bore_x)
        for index in range(1, blend_points + 1):
            angle = left_end * index / blend_points
            outline.append(
                (-centre_x + fillet * math.cos(angle), centre_y + fillet * math.sin(angle))
            )

        canvas_points: list[float] = []
        for outline_point in outline:
            canvas_points.extend(point(*outline_point))
        self.canvas.create_line(canvas_points, fill="black", width=3, joinstyle=tk.ROUND)

    def _port_status(self, state: PrescribedState, angles: CharacteristicAngles) -> str:
        """Report the port phase and the open discharge-port area.

        Only the cheap geometric quantities go here. The true-geometry
        chamber volumes need an area integration per crank angle, which is
        far too slow to run every frame.
        """

        normalized = state.crank_angle_rad % (2.0 * math.pi)
        if normalized < angles.suction_open_rad:
            phase = "suction sealed (re-expansion)"
        elif normalized < angles.compression_start_rad:
            phase = "suction"
        elif normalized < angles.discharge_open_rad:
            phase = "compression"
        elif normalized < angles.discharge_close_rad:
            phase = "discharge port open"
        else:
            phase = "recompression (port shut)"
        open_area_mm2 = port_open_area_m2(self.geometry, state.crank_angle_rad) / (MM * MM)
        text = (
            f"Port phase: {phase}    "
            f"Discharge port open area: {open_area_mm2:6.2f} mm^2    "
            f"phi {math.degrees(angles.suction_open_rad):.1f}  "
            f"beta {math.degrees(angles.compression_start_rad):.1f}  "
            f"port {math.degrees(angles.discharge_open_rad):.1f}-"
            f"{math.degrees(angles.discharge_close_rad):.1f} deg"
        )
        return text[:PORT_STATUS_WIDTH].ljust(PORT_STATUS_WIDTH)

    def _set_status(self, state: PrescribedState, angles: CharacteristicAngles) -> None:
        rotor_x, rotor_y = state.rotor_center_m
        _, cutout_y = state.cutout_center_m
        display_rpm = self.geometry.speed_rpm * self.slow_factor
        self.status_var.set(
            f"Clockwise crank angle: {math.degrees(state.crank_angle_rad):7.2f} deg    "
            f"Rotor center: ({rotor_x / MM:+.3f}, {rotor_y / MM:+.3f}) mm    "
            f"Cutout center y: {cutout_y / MM:.3f} mm    "
            f"Vane length: {state.vane_length_m / MM:.3f} mm\n"
            f"Physical speed: {self.geometry.frequency_hz:g} Hz = "
            f"{self.geometry.speed_rpm:g} rpm    "
            f"Displayed speed: {display_rpm:g} rpm    "
            f"Rotor orientation: {math.degrees(state.rotor_orientation_rad):.2f} deg\n"
            f"{self._chamber_status(state)}\n"
            f"{self._port_status(state, angles)}"
        )

    def _chamber_status(self, state: PrescribedState) -> str:
        """Report the circular-rotor chamber split without changing the drawing.

        Every branch is padded to the same character count so the fixed-width
        status line never changes the window's requested size.
        """

        try:
            areas = chamber_areas(self.geometry, state.crank_angle_rad)
            pressures = chamber_pressures(self.geometry, state.crank_angle_rad)
        except ChamberSeparationError as error:
            if in_seal_over_window(self.geometry, state.crank_angle_rad):
                half_deg = math.degrees(seal_over_half_angle_rad(self.geometry))
                mixing = chamber_pressures(self.geometry, state.crank_angle_rad)
                text = (
                    f"Seal-over mixing (+/-{half_deg:.2f} deg): both chambers at "
                    f"{mixing.discharge_pressure_pa / 1.0e6:.2f} MPa, ramping to suction"
                )
            else:
                text = f"Chamber split unavailable: {error}"
            return text[:STATUS_WIDTH].ljust(STATUS_WIDTH)
        height_m = self.geometry.cylinder_height_m
        suction_mm2 = areas.suction_area_m2 / (MM * MM)
        discharge_mm2 = areas.discharge_area_m2 / (MM * MM)
        suction_cm3 = areas.suction_area_m2 * height_m * 1.0e6
        discharge_cm3 = areas.discharge_area_m2 * height_m * 1.0e6
        valve = "valve open" if pressures.discharge_valve_open else ""
        return (
            f"Suction (IN): {suction_mm2:8.2f} mm^2 = {suction_cm3:6.3f} cm^3 "
            f"@ {pressures.suction_pressure_pa / 1.0e6:5.2f} MPa    "
            f"Discharge (OUT): {discharge_mm2:8.2f} mm^2 = {discharge_cm3:6.3f} cm^3 "
            f"@ {pressures.discharge_pressure_pa / 1.0e6:5.2f} MPa  {valve:<12s}"
        )


def main() -> int:
    """Launch the test GUI."""

    root = tk.Tk()
    RotaryCompressorApp(root)
    root.mainloop()
    return 0
