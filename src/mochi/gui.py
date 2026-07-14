"""Tkinter test GUI for prescribed rotary-compressor motion."""

from __future__ import annotations

import math
import time
import tkinter as tk
from collections.abc import Callable
from tkinter import ttk

from mochi.kinematics import MM, PrescribedState, RotaryGeometry, port_position, prescribed_state

CanvasTransform = Callable[[float, float], tuple[float, float]]


INLET_ANGLE_FROM_Y_DEG = 30.0
OUTLET_ANGLE_FROM_Y_DEG = -30.0
FRAME_DELAY_MS = 16
DEFAULT_SLOW_FACTOR = 0.01


class RotaryCompressorApp:
    """Animate the supplied rotary-compressor geometry on a Tk canvas."""

    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("mochi - rotary compressor prescribed-motion test")
        self.root.minsize(980, 700)

        self.geometry = RotaryGeometry.default()
        self.slow_factor = DEFAULT_SLOW_FACTOR
        self.crank_angle_rad = 0.0
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

        button_row = ttk.Frame(controls)
        button_row.grid(row=len(fields) + 1, column=0, columnspan=2, sticky=tk.EW, pady=(8, 4))
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
            row=len(fields) + 2,
            column=0,
            columnspan=2,
            sticky=tk.EW,
            pady=10,
        )
        ttk.Label(
            controls,
            text=(
                "Port markers are schematic:\n"
                "inlet +30 deg and outlet -30 deg\n"
                "from the global +y axis.\n\n"
                "Slow factor 0.01 displays\n"
                "motion at 1/100 real speed."
            ),
            justify=tk.LEFT,
        ).grid(row=len(fields) + 3, column=0, columnspan=2, sticky=tk.W)

        ttk.Label(
            controls,
            textvariable=self.error_var,
            foreground="#b42318",
            wraplength=265,
            justify=tk.LEFT,
        ).grid(row=len(fields) + 4, column=0, columnspan=2, sticky=tk.W, pady=(12, 0))

        view = ttk.Frame(outer)
        view.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.canvas = tk.Canvas(
            view,
            background="white",
            highlightthickness=1,
            highlightbackground="#c7ccd4",
        )
        self.canvas.pack(fill=tk.BOTH, expand=True)
        ttk.Label(
            view,
            textvariable=self.status_var,
            anchor=tk.W,
            justify=tk.LEFT,
            padding=(4, 8),
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

    def _read_inputs(self) -> tuple[RotaryGeometry, float]:
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
            )
            slow_factor = float(self.slow_factor_var.get())
        except ValueError as error:
            raise ValueError("All inputs must be valid numbers.") from error
        geometry.validate()
        if not math.isfinite(slow_factor) or not 0.0 < slow_factor <= 1.0:
            raise ValueError("Display slow factor must be greater than 0 and at most 1.")
        return geometry, slow_factor

    def _apply_inputs(self) -> None:
        try:
            geometry, slow_factor = self._read_inputs()
        except ValueError as error:
            self.error_var.set(str(error))
            return
        self.geometry = geometry
        self.slow_factor = slow_factor
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
            self.crank_angle_rad = (
                self.crank_angle_rad
                + self.geometry.angular_speed_rad_s * self.slow_factor * elapsed_s
            ) % (2.0 * math.pi)
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
        self._draw_cylinder_vane_outline(point, geometry, state.vane_tip_m)
        self._draw_center_guides(point, geometry, state)

        # Port ticks annotate the continuous cylinder outline without masking it.
        self._draw_port(point, geometry, INLET_ANGLE_FROM_Y_DEG, "IN")
        self._draw_port(point, geometry, OUTLET_ANGLE_FROM_Y_DEG, "OUT")

        self._set_status(state)

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

    def _draw_port(
        self,
        point: CanvasTransform,
        geometry: RotaryGeometry,
        angle_deg: float,
        label: str,
    ) -> None:
        inner = port_position(geometry.cylinder_radius_m - 1.0 * MM, angle_deg)
        outer = port_position(geometry.cylinder_radius_m + 2.5 * MM, angle_deg)
        label_point = port_position(geometry.cylinder_radius_m + 6.0 * MM, angle_deg)
        self.canvas.create_line(*point(*inner), *point(*outer), fill="black", width=2)
        label_px = point(*label_point)
        self.canvas.create_text(
            *label_px, text=label, fill="black", font=("TkDefaultFont", 10, "bold")
        )

    def _draw_rotor_outline(
        self,
        point: CanvasTransform,
        geometry: RotaryGeometry,
        state: PrescribedState,
    ) -> None:
        """Draw one rotor contour.

        The lips are filleted and the opening ends at the circular cutout.

        """
        rotor_x, rotor_y = state.rotor_center_m
        cutout_x, cutout_y = state.cutout_center_m
        axis_x = (cutout_x - rotor_x) / geometry.cutout_offset_m
        axis_y = (cutout_y - rotor_y) / geometry.cutout_offset_m
        normal_x = -axis_y
        normal_y = axis_x

        def world_point(axial_m: float, transverse_m: float) -> tuple[float, float]:
            return (
                rotor_x + axial_m * axis_x + transverse_m * normal_x,
                rotor_y + axial_m * axis_y + transverse_m * normal_y,
            )

        clearance = min(
            0.5 * MM,
            0.25 * (2.0 * geometry.cutout_radius_m - geometry.vane_width_m),
        )
        half_slot_width = 0.5 * geometry.vane_width_m + clearance
        rotor_radius = geometry.rotor_radius_m
        cutout_radius = geometry.cutout_radius_m
        cutout_offset = geometry.cutout_offset_m

        fillet_radius = min(
            1.5 * MM,
            0.25 * (rotor_radius - half_slot_width),
        )
        fillet_center_transverse = half_slot_width + fillet_radius
        fillet_center_radius = rotor_radius - fillet_radius
        fillet_center_axial = math.sqrt(fillet_center_radius**2 - fillet_center_transverse**2)
        fillet_angle = math.atan2(fillet_center_transverse, fillet_center_axial)
        circle_angle = math.asin(half_slot_width / cutout_radius)
        circle_axial = math.sqrt(cutout_radius**2 - half_slot_width**2)

        local_outline: list[tuple[float, float]] = []
        for index in range(121):
            angle = fillet_angle + (2.0 * math.pi - 2.0 * fillet_angle) * index / 120
            local_outline.append((rotor_radius * math.cos(angle), rotor_radius * math.sin(angle)))

        for index in range(21):
            angle = -fillet_angle + (0.5 * math.pi + fillet_angle) * index / 20
            local_outline.append(
                (
                    fillet_center_axial + fillet_radius * math.cos(angle),
                    -fillet_center_transverse + fillet_radius * math.sin(angle),
                )
            )

        local_outline.append((cutout_offset + circle_axial, -half_slot_width))
        for index in range(81):
            angle = -circle_angle + (-2.0 * math.pi + 2.0 * circle_angle) * index / 80
            local_outline.append(
                (
                    cutout_offset + cutout_radius * math.cos(angle),
                    cutout_radius * math.sin(angle),
                )
            )
        local_outline.append((fillet_center_axial, half_slot_width))

        for index in range(21):
            angle = -0.5 * math.pi + (0.5 * math.pi + fillet_angle) * index / 20
            local_outline.append(
                (
                    fillet_center_axial + fillet_radius * math.cos(angle),
                    fillet_center_transverse + fillet_radius * math.sin(angle),
                )
            )

        canvas_points: list[float] = []
        for local_point in local_outline:
            canvas_points.extend(point(*world_point(*local_point)))
        self.canvas.create_polygon(
            canvas_points,
            fill="#d9d9d9",
            outline="black",
            width=3,
            joinstyle=tk.ROUND,
        )

    def _draw_cylinder_vane_outline(
        self,
        point: CanvasTransform,
        geometry: RotaryGeometry,
        vane_tip: tuple[float, float],
    ) -> None:
        """Draw the cylinder and changing-length vane as one transparent outline."""

        radius = geometry.cylinder_radius_m
        half_vane_width = 0.5 * geometry.vane_width_m
        attachment_y = math.sqrt(radius**2 - half_vane_width**2)
        right_angle = math.atan2(attachment_y, half_vane_width)
        left_angle = math.pi - right_angle

        outline: list[tuple[float, float]] = []
        for index in range(181):
            angle = left_angle + (2.0 * math.pi + right_angle - left_angle) * index / 180
            outline.append((radius * math.cos(angle), radius * math.sin(angle)))
        outline.extend(
            (
                (half_vane_width, vane_tip[1]),
                (-half_vane_width, vane_tip[1]),
                (-half_vane_width, attachment_y),
            )
        )

        canvas_points: list[float] = []
        for outline_point in outline:
            canvas_points.extend(point(*outline_point))
        self.canvas.create_line(canvas_points, fill="black", width=3, joinstyle=tk.ROUND)

    def _set_status(self, state: PrescribedState) -> None:
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
            f"Rotor orientation: {math.degrees(state.rotor_orientation_rad):.2f} deg"
        )


def main() -> int:
    """Launch the test GUI."""

    root = tk.Tk()
    RotaryCompressorApp(root)
    root.mainloop()
    return 0
