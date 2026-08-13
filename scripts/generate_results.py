"""Regenerate the git-ignored ``results/`` figures from the current model.

This script is the single source of truth for the ``results/`` gallery. Every
figure it writes is rendered from ``mochi``'s current geometry and physics, so
running it after a geometry or model change refreshes the figures instead of
leaving a picture of a superseded shape behind.

Two figures are produced:

* ``chamber_pressures_vs_crank.png`` -- suction and compression chamber
  pressures over one revolution from the port-timed, true-geometry-volume
  rule (:func:`mochi.chambers.port_timed_pressures`, PHYSICS.md section 3.4),
  including the re-expansion and recompression phases the circular-rotor
  approximation cannot represent.
* ``rotor_motion.gif`` -- the prescribed mechanism through one revolution in
  its latest confirmed geometry: the asymmetric rotor mouth, the R2.1 mm
  vane-root fillets, and both swing-bush pieces, with the same port-timed
  pressure trace tracked alongside.

The ``results/`` directory is git-ignored and never shared, so only the latest
render of each figure is worth keeping. ``--prune`` enforces that by deleting
any image in ``results/`` that this script does not currently produce (see
AGENTS.md, "Generated and local files"). Without ``--prune`` the stale images
are only listed, not removed.

Requires the plotting extra::

    python -m pip install -e ".[viz]"
    python scripts/generate_results.py            # regenerate, list stale
    python scripts/generate_results.py --prune    # regenerate and delete stale
"""

from __future__ import annotations

import argparse
from dataclasses import replace
from math import atan2, cos, degrees, hypot, log, pi, radians, sin, sqrt
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Arc, Rectangle, Wedge

from mochi.arc_film import arc_film_force
from mochi.bearing_load import mechanism_load, shaft_work_j
from mochi.bush_film import (
    AXIAL_CLEARANCE_M,
    LUBRICANT_VISCOSITY_PA_S,
    curved_slide_velocity,
    film_state,
    film_thicknesses_m,
    flat_contact_length_m,
    flat_slide_velocity,
    friction_power_cycle_w,
)
from mochi.bush_outline import bush_piece_outline_mm as _bush_outline_mm
from mochi.bush_outline import vane_outline_mm as _vane_outline_mm
from mochi.chamber_volume import AxialBands, SwingBush, clearance_volume_m3
from mochi.chambers import (
    DISCHARGE_PORT_PRESSURE_PA,
    DISCHARGE_VALVE_RISE_FRACTION,
    SUCTION_PORT_PRESSURE_PA,
    build_cycle_trace,
    chamber_areas,
    port_timed_pressures,
    seal_over_half_angle_rad,
)
from mochi.gas_force import gas_load, gas_torque_work_j
from mochi.indicated_work import indicated_work_j
from mochi.journal_bearing import (
    JOURNAL_CLEARANCE_M,
    JOURNAL_LENGTH_M,
    JOURNAL_RADIUS_M,
    journal_relative_speed_rad_s,
)
from mochi.kinematics import (
    MM,
    RotaryGeometry,
    port_position,
    prescribed_state,
    vane_fillet_geometry,
)
from mochi.leakage import SUCTION_DENSITY_KG_M3, leaky_cycle
from mochi.long_bearing import long_bearing_load
from mochi.ocvirk_bearing import eccentricity_cycle, short_bearing_force
from mochi.ports import (
    characteristic_angles,
    discharge_window,
    port_open_area_m2,
    suction_window,
)
from mochi.reed_valve import valved_cycle
from mochi.reynolds_1d import solve_short_bearing_1d
from mochi.rotor_bush_dynamics import integrate_rotor_bush_orbit
from mochi.rotor_cylinder import contact_normal_force_n, rotor_cylinder_friction_power_w
from mochi.rotor_dynamics import integrate_rotor_orbit
from mochi.rotor_profile import rotor_contour
from mochi.slider_film import flat_slider_film
from mochi.tecplot import point_zone, write_dat
from mochi.thermo_check import isentropic_cross_check
from mochi.true_gas_force import true_gas_load, true_gas_torque_work_j

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
# Raw model data (Tecplot ASCII .dat) for research post-processing. The PNGs under
# results/ are for illustration; the numbers behind them live here (see DATASETS).
DATA_DIR = RESULTS_DIR / "data"

MPA = 1.0e6
PA_TO_MPA = 1.0 / MPA

# Fill colours, matched to the earlier hand-made figures.
GAS_COLOR = "#e8934f"
ROTOR_COLOR = "#d9d9d9"
BUSH_COLOR = "#b9c0cb"
VANE_COLOR = "#5a5a5a"
GUIDE_COLOR = "#5c6068"
SUCTION_COLOR = "#3d7fc1"
COMPRESSION_COLOR = "#e8752a"


def _use_korean_font() -> None:
    """Pick an installed CJK font so the Korean labels render as glyphs."""

    available = {f.name for f in font_manager.fontManager.ttflist}
    for name in ("Malgun Gothic", "Gulim", "Dotum", "Batang", "NanumGothic"):
        if name in available:
            plt.rcParams["font.family"] = name
            break
    plt.rcParams["axes.unicode_minus"] = False


def _save(fig, path: Path) -> None:
    """Create the parent folder, tight-layout, save, and close the figure."""

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _legend_above(ax, ncol: int, fontsize: int = 9) -> None:
    """Place the legend above the axes (outside the plot) so it never covers data."""

    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.0),
        ncol=ncol,
        fontsize=fontsize,
        framealpha=0.9,
    )


def _reed_note(ax, x: float = 0.02, y: float = 0.98, ha: str = "left", va: str = "top") -> None:
    """Small grey note: these mechanical loads use the baseline (reed-valve-free)
    basis — the isentropic-ideal indicated work. The discharge reed-valve
    overpressure (Section 3.8) is a separate performance term, not propagated."""

    ax.text(
        x,
        y,
        "리드밸브 미포함 (이상 지시 기준)",
        transform=ax.transAxes,
        ha=ha,
        va=va,
        fontsize=8,
        color="#888",
        zorder=20,
    )


# --------------------------------------------------------------------------
# Shared pressure trace
# --------------------------------------------------------------------------


class PressureCurve:
    """Port-timed chamber pressures sampled over one revolution."""

    def __init__(self, geometry: RotaryGeometry, samples: int = 721) -> None:
        trace = build_cycle_trace(geometry)
        self.geometry = geometry
        self.trace = trace  # reused by the P-V and bush-film figures
        self.angles_deg: list[float] = []
        self.suction_mpa: list[float] = []
        self.compression_mpa: list[float] = []
        phases: list[str] = []
        for index in range(samples):
            angle_rad = 2.0 * pi * index / (samples - 1)
            result = port_timed_pressures(geometry, angle_rad, trace=trace)
            self.angles_deg.append(degrees(angle_rad))
            self.suction_mpa.append(result.suction_pressure_pa * PA_TO_MPA)
            self.compression_mpa.append(result.compression_pressure_pa * PA_TO_MPA)
            phases.append(result.phase)

        angles = characteristic_angles(geometry)
        self.phi_deg = degrees(angles.suction_open_rad)
        self.beta_deg = degrees(angles.compression_start_rad)
        self.discharge_close_deg = degrees(angles.discharge_close_rad)
        self.seal_over_entry_deg = 360.0 - degrees(seal_over_half_angle_rad(geometry))
        # The valve angle is a pressure condition, not a constant, so read it
        # back from where the compression phase first becomes delivery.
        self.valve_open_deg = next(
            (self.angles_deg[i] for i, phase in enumerate(phases) if phase == "discharge"),
            self.beta_deg,
        )
        self.residual_mpa = (
            port_timed_pressures(geometry, 0.0, trace=trace).residual_pressure_pa * PA_TO_MPA
        )
        self.opening_mpa = (
            DISCHARGE_PORT_PRESSURE_PA * (1.0 + DISCHARGE_VALVE_RISE_FRACTION) * PA_TO_MPA
        )
        self.suction_port_mpa = SUCTION_PORT_PRESSURE_PA * PA_TO_MPA
        self.discharge_port_mpa = DISCHARGE_PORT_PRESSURE_PA * PA_TO_MPA
        self.clearance_cm3 = clearance_volume_m3(geometry) * 1.0e6

        # Reed-valve overpressure (Section 3.8) — a **separate performance term**,
        # kept only for the P-V diagram's labelled overlay; NOT propagated into the
        # baseline loads (force/torque/bearing/shaft all use the ideal indicated).
        valved = valved_cycle(geometry, trace=trace, samples=2880)
        self.baseline_indicated_power_w = valved.baseline_indicated_work_j * geometry.frequency_hz
        self.overpressure_power_w = valved.overpressure_power_w
        self.valve_indicated_power_w = valved.valve_indicated_power_w

    def value_at(self, angle_deg: float) -> tuple[float, float]:
        """Nearest sampled (suction, compression) pressure for a crank angle."""

        return self.suction_mpa[self._index(angle_deg)], self.compression_mpa[
            self._index(angle_deg)
        ]

    def _index(self, angle_deg: float) -> int:
        span = self.angles_deg[-1] - self.angles_deg[0]
        index = round((angle_deg % 360.0 - self.angles_deg[0]) / span * (len(self.angles_deg) - 1))
        return min(max(index, 0), len(self.angles_deg) - 1)


# --------------------------------------------------------------------------
# Figure 1: chamber pressures vs crank angle
# --------------------------------------------------------------------------


def render_chamber_pressures(curve: PressureCurve, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10.5, 6.0), dpi=140)
    ax.plot(
        curve.angles_deg,
        curve.suction_mpa,
        color=SUCTION_COLOR,
        linewidth=2.2,
        label="흡입 챔버 (IN 쪽)",
    )
    ax.plot(
        curve.angles_deg,
        curve.compression_mpa,
        color=COMPRESSION_COLOR,
        linewidth=2.2,
        label="압축/토출 챔버 (OUT 쪽)",
    )

    for pressure, text, valign in (
        (curve.suction_port_mpa, f"흡입 포트압 {curve.suction_port_mpa:.2f}", "bottom"),
        (curve.discharge_port_mpa, f"토출 포트압 {curve.discharge_port_mpa:.2f}", "top"),
        (curve.opening_mpa, f"밸브 개방압 {curve.opening_mpa:.2f}", "bottom"),
    ):
        ax.axhline(pressure, color="#9aa0a8", linestyle=":", linewidth=1.0)
        ax.text(357, pressure, text, va=valign, ha="right", fontsize=8, color="#5c6068")

    phase_marks = (
        (curve.phi_deg, rf"$\varphi$ 흡입 개방 {curve.phi_deg:.1f}°", "#5c6068"),
        (curve.beta_deg, f"β 압축 시작 {curve.beta_deg:.1f}°", "#c0532a"),
        (curve.valve_open_deg, f"밸브 개방 {curve.valve_open_deg:.0f}°", "#c0532a"),
        (curve.discharge_close_deg, f"δ 토출 폐쇄 {curve.discharge_close_deg:.1f}°", "#7a4fb0"),
    )
    for angle_deg, text, color in phase_marks:
        ax.axvline(angle_deg, color=color, linestyle="--", linewidth=0.9, alpha=0.55)
        ax.text(
            angle_deg + 1.5,
            0.2,
            text,
            rotation=90,
            va="bottom",
            ha="left",
            fontsize=7.5,
            color=color,
        )

    ax.annotate(
        f"재압축 피크 약 {curve.residual_mpa:.1f} MPa\n(누설 무시 -> 상한값; §3.7에서 억제)",
        xy=(curve.seal_over_entry_deg, curve.residual_mpa),
        xytext=(250, curve.residual_mpa - 0.5),
        fontsize=9,
        color="#7a4fb0",
        arrowprops={"arrowstyle": "->", "color": "#7a4fb0", "lw": 1.0},
    )
    ax.annotate(
        r"$\varphi$ 이전: 잔류 가스 재팽창",
        xy=(curve.phi_deg * 0.4, curve.residual_mpa * 0.65),
        xytext=(40, curve.residual_mpa * 0.75),
        fontsize=9,
        color=SUCTION_COLOR,
        arrowprops={"arrowstyle": "->", "color": SUCTION_COLOR, "lw": 1.0},
    )
    ax.text(
        curve.beta_deg + 8,
        1.35,
        "압축\np·V^n = 일정",
        fontsize=9,
        color="#c0532a",
        rotation=32,
    )
    ax.text(
        (curve.valve_open_deg + curve.discharge_close_deg) / 2,
        curve.discharge_port_mpa + 0.18,
        "토출(포트압 유지)",
        fontsize=9,
        color="#c0532a",
        ha="center",
    )
    ax.text(
        curve.discharge_close_deg - 2,
        6.0,
        "δ 이후: 재압축\n(포트 폐쇄)",
        fontsize=9,
        color="#7a4fb0",
        ha="right",
    )

    ax.set_xlim(0, 360)
    ax.set_ylim(0, max(11.0, curve.residual_mpa + 1.0))
    ax.set_xticks(range(0, 361, 45))
    ax.set_xlabel("크랭크 각 θ (deg, 상단 0°에서 시계방향)")
    ax.set_ylabel("챔버 압력 (MPa, 절대압)")
    ax.set_title(
        "크랭크 각별 챔버 압력 — R410A, n = 1.07, 밸브 과압 5%\n"
        f"포트-타이밍 실형상 체적 규칙 (PHYSICS.md §3.4, 클리어런스 {curve.clearance_cm3:.3f} cm³)",
        color="none",
    )
    ax.grid(True, color="#e2e5ea", linewidth=0.8)
    _legend_above(ax, ncol=2)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# --------------------------------------------------------------------------
# Figure 2: rotor motion animation
# --------------------------------------------------------------------------


def _phase_label(geometry: RotaryGeometry, angle_deg: float) -> str:
    angles = characteristic_angles(geometry)
    half_deg = degrees(seal_over_half_angle_rad(geometry))
    normalized = angle_deg % 360.0
    if normalized < half_deg or normalized > 360.0 - half_deg:
        return "밀봉 전환(혼합)"
    if normalized < degrees(angles.suction_open_rad):
        return "흡입 밀봉(재팽창)"
    if normalized < degrees(angles.compression_start_rad):
        return "흡입"
    if normalized < degrees(angles.discharge_open_rad):
        return "압축"
    if normalized < degrees(angles.discharge_close_rad):
        return "토출 포트 개방"
    return "재압축(포트 폐쇄)"


def _draw_port_arc(
    ax,
    geometry: RotaryGeometry,
    start_deg: float,
    end_deg: float,
    label: str,
    *,
    color: str = "black",
    label_color: str | None = None,
    label_offset_mm: float = 7.0,
    tick_labels: tuple[tuple[float, str], ...] | None = None,
) -> None:
    """Stroke one port as an arc on the bore, with end ticks and a label.

    ``tick_labels`` adds a leader and a degree readout at each named angle,
    used by the dimensioned ``port_geometry`` figure; the animation leaves it
    ``None`` and only draws the IN/OUT arcs.
    """

    radius = geometry.cylinder_radius_m + 2.0 * MM
    span_deg = (end_deg - start_deg) % 360.0
    steps = max(int(span_deg), 2) + 1
    xs: list[float] = []
    ys: list[float] = []
    for index in range(steps):
        angle_deg = start_deg + span_deg * index / (steps - 1)
        x_m, y_m = port_position(radius, angle_deg)
        xs.append(x_m / MM)
        ys.append(y_m / MM)
    ax.plot(xs, ys, color=color, linewidth=3.0, solid_capstyle="round")
    for angle_deg in (start_deg, end_deg):
        inner = port_position(geometry.cylinder_radius_m - 1.0 * MM, angle_deg)
        outer = port_position(radius + 1.5 * MM, angle_deg)
        ax.plot(
            [inner[0] / MM, outer[0] / MM],
            [inner[1] / MM, outer[1] / MM],
            color=color,
            linewidth=1.6,
        )
    for angle_deg, text in tick_labels or ():
        tick = port_position(radius + 1.5 * MM, angle_deg)
        seat = port_position(radius + 5.0 * MM, angle_deg)
        ax.plot(
            [tick[0] / MM, seat[0] / MM], [tick[1] / MM, seat[1] / MM], color=color, linewidth=1.2
        )
        label_pos = port_position(radius + 6.5 * MM, angle_deg)
        ax.text(
            label_pos[0] / MM,
            label_pos[1] / MM,
            text,
            ha="center",
            va="center",
            fontweight="bold",
            fontsize=10,
            color=color,
        )
    if label:
        label_x, label_y = port_position(
            geometry.cylinder_radius_m + label_offset_mm * MM, start_deg + 0.5 * span_deg
        )
        ax.text(
            label_x / MM,
            label_y / MM,
            label,
            ha="center",
            va="center",
            fontweight="bold",
            fontsize=11,
            color=label_color or color,
        )


def render_rotor_motion(curve: PressureCurve, path: Path, frames: int = 120) -> None:
    geometry = curve.geometry
    bore_mm = geometry.cylinder_radius_m / MM
    limit = bore_mm + 9.0

    fig = plt.figure(figsize=(11.2, 5.4), dpi=110)
    ax_left = fig.add_axes([0.02, 0.06, 0.46, 0.82])
    ax_right = fig.add_axes([0.57, 0.14, 0.40, 0.72])
    frame_title = fig.text(0.25, 0.90, "", ha="center", fontsize=11)

    # Static right panel: the baseline port-timed pressure trace, with a marker.
    ax_right.plot(
        curve.angles_deg, curve.suction_mpa, color=SUCTION_COLOR, linewidth=1.8, label="흡입 챔버"
    )
    ax_right.plot(
        curve.angles_deg,
        curve.compression_mpa,
        color=COMPRESSION_COLOR,
        linewidth=1.8,
        label="압축/토출 챔버 (이상, 리드밸브 미포함)",
    )
    ax_right.set_xlim(0, 360)
    ax_right.set_ylim(0, max(11.0, curve.residual_mpa + 1.0))
    ax_right.set_xticks(range(0, 361, 60))
    ax_right.set_xlabel("크랭크 각 θ (deg)")
    ax_right.set_ylabel("챔버 압력 (MPa, 절대압)")
    ax_right.grid(True, color="#e2e5ea", linewidth=0.7)
    # Legend above the axes (outside the plot) so it never covers the traces.
    ax_right.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.0),
        ncol=2,
        fontsize=8,
        framealpha=0.9,
        borderaxespad=0.2,
    )

    # Shade the seal-over window that the indicated-work / gas-torque integrals
    # skip. Integration stops at the seal-over entry (2*pi - alpha ~ 354 deg) --
    # exactly where the recompression pressure peaks -- so the merged region past
    # it, where the true-geometry chamber split degenerates, is excluded.
    seal_half_deg = 360.0 - curve.seal_over_entry_deg
    for span_start, span_end in ((0.0, seal_half_deg), (curve.seal_over_entry_deg, 360.0)):
        ax_right.axvspan(span_start, span_end, color=MERGED_COLOR, alpha=0.22, zorder=0)
    ax_right.axvline(
        curve.seal_over_entry_deg,
        color=MERGED_COLOR,
        linestyle=(0, (4, 3)),
        linewidth=1.3,
        zorder=2,
    )
    ax_right.text(
        0.985,
        0.70,
        f"seal-over ±{seal_half_deg:.2f}°\n"
        f"θ≥{curve.seal_over_entry_deg:.0f}° 적분 제외\n(챔버 병합)",
        transform=ax_right.transAxes,
        ha="right",
        va="top",
        fontsize=7.5,
        color="#5a3f74",
        bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": MERGED_COLOR, "alpha": 0.9},
    )

    marker_line = ax_right.axvline(0.0, color=GUIDE_COLOR, linewidth=1.0)
    (suction_dot,) = ax_right.plot([], [], "o", color=SUCTION_COLOR, markersize=6)
    (compression_dot,) = ax_right.plot([], [], "o", color=COMPRESSION_COLOR, markersize=6)

    def draw_left(angle_deg: float) -> None:
        ax_left.clear()
        ax_left.set_aspect("equal")
        ax_left.set_xlim(-limit, limit)
        ax_left.set_ylim(-limit, limit)
        ax_left.axis("off")
        angle_rad = radians(angle_deg)
        state = prescribed_state(geometry, angle_rad)

        # Gas fills the whole bore; the solid parts are painted back over it.
        # The rotor-cylinder contact (at bore angle = crank angle) and the vane
        # split the crescent into the suction chamber (IN, blue) and the
        # discharge chamber (OUT, orange); through the seal-over window near top
        # dead center the two merge into one region (purple). Bore angle runs
        # from +y clockwise, so matplotlib's CCW-from-+x wedge angle is 90 - it.
        half_deg = degrees(seal_over_half_angle_rad(geometry))
        if angle_deg < half_deg or angle_deg > 360.0 - half_deg:
            ax_left.add_patch(
                plt.Circle((0.0, 0.0), bore_mm, facecolor=MERGED_COLOR, edgecolor="none", zorder=1)
            )
        else:
            ax_left.add_patch(
                Wedge(
                    (0.0, 0.0), bore_mm, 90.0 - angle_deg, 90.0, facecolor=SUCTION_COLOR, zorder=1
                )
            )
            ax_left.add_patch(
                Wedge(
                    (0.0, 0.0),
                    bore_mm,
                    -270.0,
                    90.0 - angle_deg,
                    facecolor=COMPRESSION_COLOR,
                    zorder=1,
                )
            )

        contour = rotor_contour(geometry, angle_rad)
        material_x = [p[0] / MM for p in contour.material]
        material_y = [p[1] / MM for p in contour.material]
        ax_left.fill(material_x, material_y, facecolor=ROTOR_COLOR, edgecolor="none", zorder=2)

        groove_x, groove_y = state.cutout_center_m
        bush = SwingBush()
        for side in (1.0, -1.0):
            centre_x = groove_x + side * bush.piece_shift_m
            xs, ys = _bush_outline_mm(centre_x, groove_y, side)
            ax_left.fill(xs, ys, facecolor=BUSH_COLOR, edgecolor="black", linewidth=0.8, zorder=3)

        vane_x, vane_y = _vane_outline_mm(geometry, state.vane_tip_m[1])
        ax_left.fill(
            vane_x, vane_y, facecolor=VANE_COLOR, edgecolor="black", linewidth=1.2, zorder=4
        )

        # Rotor-cylinder contact point T (the moving IN/OUT chamber divider).
        ax_left.plot(
            bore_mm * sin(angle_rad),
            bore_mm * cos(angle_rad),
            "o",
            color="#c0392b",
            markersize=6,
            zorder=6,
        )

        # Bore outline and the real rotor material edges.
        ax_left.add_patch(
            plt.Circle(
                (0.0, 0.0), bore_mm, facecolor="none", edgecolor="black", linewidth=2.0, zorder=5
            )
        )
        for edge in (contour.od_arc, contour.inlet_flat, contour.mouth_path):
            ax_left.plot(
                [p[0] / MM for p in edge],
                [p[1] / MM for p in edge],
                color="black",
                linewidth=2.2,
                zorder=5,
            )

        # Rotor-centre orbit and the C / R markers.
        eccentricity_mm = geometry.eccentricity_m / MM
        ax_left.add_patch(
            plt.Circle(
                (0.0, 0.0),
                eccentricity_mm,
                facecolor="none",
                edgecolor=GUIDE_COLOR,
                linewidth=1.0,
                linestyle=(0, (5, 4)),
                zorder=6,
            )
        )
        for (cx_m, cy_m), text in (((0.0, 0.0), "C"), (state.rotor_center_m, "R")):
            cx, cy = cx_m / MM, cy_m / MM
            ax_left.plot(
                cx, cy, marker="+", color=GUIDE_COLOR, markersize=9, markeredgewidth=2, zorder=7
            )
            ax_left.text(
                cx + 1.4, cy + 1.4, text, color=GUIDE_COLOR, fontweight="bold", fontsize=9, zorder=7
            )

        suction = suction_window(geometry)
        discharge = discharge_window(geometry)
        _draw_port_arc(
            ax_left, geometry, degrees(suction.start_rad), degrees(suction.end_rad), "IN"
        )
        _draw_port_arc(
            ax_left, geometry, degrees(discharge.start_rad), degrees(discharge.end_rad), "OUT"
        )

    def update(frame: int):
        angle_deg = 360.0 * frame / frames
        draw_left(angle_deg)
        suction_mpa, compression_mpa = curve.value_at(angle_deg)
        marker_line.set_xdata([angle_deg, angle_deg])
        suction_dot.set_data([angle_deg], [suction_mpa])
        compression_dot.set_data([angle_deg], [compression_mpa])
        in_seal_over = angle_deg < seal_half_deg or angle_deg > curve.seal_over_entry_deg
        seal_note = "   ·   seal-over 병합 (적분 제외)" if in_seal_over else ""
        frame_title.set_text(
            f"θ = {angle_deg:5.1f}°   ·   {_phase_label(geometry, angle_deg)}   ·   "
            f"IN {suction_mpa:.2f} / OUT {compression_mpa:.2f} MPa{seal_note}"
        )
        return marker_line, suction_dot, compression_dot, frame_title

    animation = FuncAnimation(fig, update, frames=frames, blit=False)
    animation.save(path, writer=PillowWriter(fps=20))
    plt.close(fig)


# --------------------------------------------------------------------------
# Figure 3: dimensioned port geometry at top dead center
# --------------------------------------------------------------------------

DISCHARGE_PORT_COLOR = "#c0392b"
SUCTION_PORT_COLOR = "#1a9979"


def render_port_geometry(geometry: RotaryGeometry, path: Path) -> None:
    """Final rotor, swing bush, and dimensioned port windows at theta = 0."""

    bore_mm = geometry.cylinder_radius_m / MM
    state = prescribed_state(geometry, 0.0)
    angles = characteristic_angles(geometry)
    suction = suction_window(geometry)
    discharge = discharge_window(geometry)

    fig, ax = plt.subplots(figsize=(9.0, 9.4), dpi=130)
    ax.set_aspect("equal")
    ax.set_xlim(-bore_mm - 8.0, bore_mm + 8.0)
    ax.set_ylim(-bore_mm - 10.0, bore_mm + 22.0)
    ax.axis("off")

    fig.suptitle(
        f"Final rotor, swing bush, and port windows  (crank angle 0°)\n"
        f"Suction φ={degrees(angles.suction_open_rad):.1f}°→"
        f"β={degrees(angles.compression_start_rad):.1f}°   |   "
        f"Discharge {degrees(angles.discharge_open_rad):.1f}°→"
        f"{degrees(angles.discharge_close_rad):.1f}° "
        f"(γ={geometry.discharge_port_span_deg:.1f}°, "
        f"shut δ={geometry.recompression_angle_deg:.1f}° before TDC)",
        fontsize=12,
    )

    # Rotor material, swing bush, and the white vane slot (no gas fill here:
    # this is a geometry drawing, not a cycle frame).
    contour = rotor_contour(geometry, 0.0)
    ax.fill(
        [p[0] / MM for p in contour.material],
        [p[1] / MM for p in contour.material],
        facecolor=ROTOR_COLOR,
        edgecolor="none",
        zorder=2,
    )
    groove_x, groove_y = state.cutout_center_m
    bush = SwingBush()
    for side in (1.0, -1.0):
        xs, ys = _bush_outline_mm(groove_x + side * bush.piece_shift_m, groove_y, side)
        ax.fill(xs, ys, facecolor=BUSH_COLOR, edgecolor="black", linewidth=0.9, zorder=3)
    vane_x, vane_y = _vane_outline_mm(geometry, state.vane_tip_m[1])
    ax.fill(vane_x, vane_y, facecolor="white", edgecolor="black", linewidth=1.4, zorder=4)

    ax.add_patch(
        plt.Circle(
            (0.0, 0.0), bore_mm, facecolor="none", edgecolor="black", linewidth=2.2, zorder=5
        )
    )
    for edge in (contour.od_arc, contour.inlet_flat, contour.mouth_path):
        ax.plot(
            [p[0] / MM for p in edge],
            [p[1] / MM for p in edge],
            color="black",
            linewidth=2.2,
            zorder=5,
        )

    for (cx_m, cy_m), text in (((0.0, 0.0), "C"), (state.rotor_center_m, "R")):
        cx, cy = cx_m / MM, cy_m / MM
        ax.plot(cx, cy, marker="+", color=GUIDE_COLOR, markersize=10, markeredgewidth=2, zorder=7)
        ax.text(
            cx + 1.6, cy + 1.0, text, color=GUIDE_COLOR, fontweight="bold", fontsize=11, zorder=7
        )

    # Top dead center reference line and the clockwise-angle sense arrow.
    ax.plot(
        [0.0, 0.0],
        [bore_mm + 1.0, bore_mm + 13.0],
        color="#7a7f88",
        linestyle=(0, (5, 4)),
        linewidth=1.2,
    )
    ax.text(
        0.0, bore_mm + 13.5, "TDC (+y, θ=0)", color="#7a7f88", fontsize=10, ha="center", va="bottom"
    )
    ax.annotate(
        "θ (cw)",
        xy=(bore_mm * 0.66, bore_mm * 0.60),
        xytext=(bore_mm + 1.0, bore_mm * 0.86),
        fontsize=10,
        color="#7a7f88",
        ha="left",
        arrowprops={"arrowstyle": "->", "color": "#7a7f88", "lw": 1.2},
    )

    # Port arcs with angle ticks; the bold headers are placed by hand so they
    # clear the TDC line and the sense arrow.
    _draw_port_arc(
        ax,
        geometry,
        degrees(discharge.start_rad),
        degrees(discharge.end_rad),
        "",
        color=DISCHARGE_PORT_COLOR,
        tick_labels=(
            (degrees(discharge.start_rad), f"{degrees(discharge.start_rad):.1f}°"),
            (degrees(discharge.end_rad), f"{degrees(discharge.end_rad):.1f}°"),
        ),
    )
    _draw_port_arc(
        ax,
        geometry,
        degrees(suction.start_rad),
        degrees(suction.end_rad),
        "",
        color=SUCTION_PORT_COLOR,
        tick_labels=(
            (degrees(suction.start_rad), f"{degrees(suction.start_rad):.1f}°"),
            (degrees(suction.end_rad), f"{degrees(suction.end_rad):.1f}°"),
        ),
    )
    ax.text(
        -bore_mm * 0.52,
        bore_mm + 9.0,
        "Discharge port\nγ span",
        color=DISCHARGE_PORT_COLOR,
        fontsize=11,
        fontweight="bold",
        ha="center",
        va="bottom",
    )
    ax.text(
        bore_mm * 0.52,
        bore_mm + 9.0,
        "Suction port\nφ→β",
        color=SUCTION_PORT_COLOR,
        fontsize=11,
        fontweight="bold",
        ha="center",
        va="bottom",
    )

    # Point at the right vane-root blend arc from the lower right.
    _, _, (bore_tangent_x, bore_tangent_y) = vane_fillet_geometry(geometry)
    ax.annotate(
        "R2.1 vane-root blend",
        xy=(bore_tangent_x / MM, bore_tangent_y / MM),
        xytext=(bore_mm * 0.40, bore_mm * 0.30),
        fontsize=10,
        color="#5c6068",
        ha="left",
        arrowprops={"arrowstyle": "->", "color": "#5c6068", "lw": 1.2},
    )

    fig.savefig(path)
    plt.close(fig)


# --------------------------------------------------------------------------
# Figure 4: top / bottom dead centre chamber definition
# --------------------------------------------------------------------------

MERGED_COLOR = "#9b72b0"


def _fill_crescent_sector(
    ax,
    geometry: RotaryGeometry,
    theta_rad: float,
    a0_deg: float,
    a1_deg: float,
    color: str,
) -> None:
    """Fill the crescent between bore and circular rotor over a bore-angle band.

    The inner edge is the circular-rotor disc (the section 3.2 approximation),
    solved along each ray from the bore centre; the outer edge is the bore.
    """

    bore = geometry.cylinder_radius_m / MM
    rotor = geometry.rotor_radius_m / MM
    ecc = geometry.eccentricity_m / MM
    steps = max(int(abs(a1_deg - a0_deg)), 2) + 1
    outer: list[tuple[float, float]] = []
    inner: list[tuple[float, float]] = []
    for index in range(steps):
        a = radians(a0_deg + (a1_deg - a0_deg) * index / (steps - 1))
        outer.append((bore * sin(a), bore * cos(a)))
        offset = a - theta_rad
        rho = ecc * cos(offset) + sqrt(max(rotor * rotor - ecc * ecc * sin(offset) ** 2, 0.0))
        inner.append((rho * sin(a), rho * cos(a)))
    poly = outer + inner[::-1]
    ax.fill([p[0] for p in poly], [p[1] for p in poly], facecolor=color, edgecolor="none", zorder=1)


def _draw_dead_centre_panel(
    ax, geometry: RotaryGeometry, theta_deg: float, title: str, merged: bool
) -> None:
    bore = geometry.cylinder_radius_m / MM
    rotor = geometry.rotor_radius_m / MM
    ecc = geometry.eccentricity_m / MM
    theta = radians(theta_deg)
    half_deg = degrees(seal_over_half_angle_rad(geometry))

    ax.set_aspect("equal")
    ax.set_xlim(-bore - 6, bore + 6)
    ax.set_ylim(-bore - 10, bore + 8)
    ax.axis("off")
    ax.set_title(title, fontsize=13)

    if merged:
        _fill_crescent_sector(ax, geometry, theta, half_deg, 360.0 - half_deg, MERGED_COLOR)
    else:
        _fill_crescent_sector(ax, geometry, theta, half_deg, theta_deg, SUCTION_COLOR)
        _fill_crescent_sector(ax, geometry, theta, theta_deg, 360.0 - half_deg, COMPRESSION_COLOR)

    orbit_x, orbit_y = ecc * sin(theta), ecc * cos(theta)
    ax.add_patch(
        plt.Circle(
            (orbit_x, orbit_y),
            rotor,
            facecolor=ROTOR_COLOR,
            edgecolor="black",
            linewidth=1.5,
            zorder=2,
        )
    )
    ax.add_patch(
        plt.Circle((0.0, 0.0), bore, facecolor="none", edgecolor="black", linewidth=2.2, zorder=3)
    )

    state = prescribed_state(geometry, theta)
    half_width = 0.5 * geometry.vane_width_m / MM
    vane_tip_y = state.vane_tip_m[1] / MM
    ax.add_patch(
        Rectangle(
            (-half_width, vane_tip_y),
            2.0 * half_width,
            bore - vane_tip_y,
            facecolor=VANE_COLOR,
            edgecolor="black",
            linewidth=1.0,
            zorder=4,
        )
    )

    # Rotor-cylinder contact point T rides the bore at the crank angle.
    ax.plot(bore * sin(theta), bore * cos(theta), "o", color="#c0392b", markersize=9, zorder=6)

    for (cx, cy), text in (((0.0, 0.0), "C"), ((orbit_x, orbit_y), "R")):
        ax.plot(cx, cy, marker="+", color=GUIDE_COLOR, markersize=10, markeredgewidth=2, zorder=5)
        ax.text(
            cx + 1.6, cy, text, color=GUIDE_COLOR, fontsize=11, ha="left", va="center", zorder=5
        )

    for angle_deg, text, ha in ((-40.0, "OUT", "right"), (40.0, "IN", "left")):
        tick_in = port_position(geometry.cylinder_radius_m - 1.0 * MM, angle_deg)
        tick_out = port_position(geometry.cylinder_radius_m + 2.5 * MM, angle_deg)
        ax.plot(
            [tick_in[0] / MM, tick_out[0] / MM],
            [tick_in[1] / MM, tick_out[1] / MM],
            color="black",
            linewidth=2,
        )
        label = port_position(geometry.cylinder_radius_m + 5.0 * MM, angle_deg)
        ax.text(
            label[0] / MM, label[1] / MM, text, fontweight="bold", fontsize=11, ha=ha, va="center"
        )


def render_tdc_bdc_definition(geometry: RotaryGeometry, path: Path) -> None:
    """Top/bottom dead-centre chamber definition (circular-rotor, section 3.2)."""

    areas = chamber_areas(geometry, pi)
    area_mm2 = areas.suction_area_m2 / (MM * MM)
    volume_cm3 = areas.suction_area_m2 * geometry.cylinder_height_m * 1.0e6

    fig, (ax_tdc, ax_bdc) = plt.subplots(1, 2, figsize=(13.0, 7.4), dpi=120)
    fig.suptitle(
        "상사점·하사점 정의 (기하 기준) — 크랭크각 θ: "
        "상단 0°에서 시계방향, 접촉점 T가 θ를 따라 회전",
        fontsize=13,
    )

    _draw_dead_centre_panel(ax_tdc, geometry, 0.0, "상사점 TDC — θ = 0°", merged=True)
    ax_tdc.annotate(
        "접촉점 T = 베인 위치\n(로터가 베인 쪽 벽면에 최근접)",
        xy=(0.0, geometry.cylinder_radius_m / MM),
        xytext=(-geometry.cylinder_radius_m / MM - 4, geometry.cylinder_radius_m / MM + 4),
        fontsize=9,
        color="#c0392b",
        ha="left",
        arrowprops={"arrowstyle": "->", "color": "#c0392b", "lw": 1.0},
    )
    ax_tdc.text(
        0.5,
        -0.02,
        "토출 부피 = 0 (토출 종료)\n흡입 챔버 밀봉 = 최대 부피 도달 → 압축 시작\n"
        f"혼합 구간 ±{degrees(seal_over_half_angle_rad(geometry)):.2f}°의 중심",
        transform=ax_tdc.transAxes,
        ha="center",
        va="top",
        fontsize=9,
    )

    _draw_dead_centre_panel(ax_bdc, geometry, 180.0, "하사점 BDC — θ = 180°", merged=False)
    ax_bdc.annotate(
        "접촉점 T = 베인 정반대\n(로터가 베인 쪽 벽면에서 최원)",
        xy=(0.0, -geometry.cylinder_radius_m / MM),
        xytext=(geometry.cylinder_radius_m / MM * 0.35, -geometry.cylinder_radius_m / MM - 3),
        fontsize=9,
        color="#c0392b",
        ha="left",
        arrowprops={"arrowstyle": "->", "color": "#c0392b", "lw": 1.0},
    )
    ax_bdc.text(
        0.5,
        -0.02,
        f"흡입(파랑) = 압축(주황) = 각 {area_mm2:.0f} mm² ({volume_cm3:.1f} cm³)\n"
        "부피 변화율 |dV/dθ| 최대인 지점 (열역학적 사점 아님)",
        transform=ax_bdc.transAxes,
        ha="center",
        va="top",
        fontsize=9,
    )

    fig.tight_layout(rect=(0.0, 0.12, 1.0, 0.95))
    fig.savefig(path)
    plt.close(fig)


# --------------------------------------------------------------------------
# Figure 5: rotor-mouth lip detail (current asymmetric R1.5 / R1.0 lips)
# --------------------------------------------------------------------------

# The default rotor_contour sampling: the inlet lip is [start] + 41-point R1.5
# arc + [joint] + 31-point R1.0 blend = 74 points at the head of mouth_path,
# followed by an 81-point groove arc, then the reversed outlet lip.
_INLET_LIP_LEN = 1 + 41 + 1 + 31
_GROOVE_LEN = 81


def render_rotor_mouth_lip_detail(geometry: RotaryGeometry, path: Path) -> None:
    """Inlet-lip G1 chain (OD flat -> R1.5 arc -> straight -> R1.0 groove blend).

    Reads the lip straight off the current asymmetric contour, so the radii and
    sweeps are whatever ``rotor_profile`` currently builds (R1.5 / R1.0), not
    the superseded R1.4 / R0.9 the earlier hand-made figure drew.
    """

    bore_mm = geometry.cylinder_radius_m / MM
    state = prescribed_state(geometry, 0.0)
    contour = rotor_contour(geometry, 0.0)
    mouth = [(x / MM, y / MM) for x, y in contour.mouth_path]
    groove_x, groove_y = (c / MM for c in state.cutout_center_m)

    inlet_touch = mouth[_INLET_LIP_LEN]  # groove touch on the inlet lip
    outlet_touch = mouth[_INLET_LIP_LEN + _GROOVE_LEN - 1]
    inlet_angle = degrees(atan2(inlet_touch[0] - groove_x, inlet_touch[1] - groove_y))

    fig = plt.figure(figsize=(12.0, 9.2), dpi=120)
    ax = fig.add_axes([0.03, 0.04, 0.60, 0.88])
    ax.set_aspect("equal")
    ax.set_xlim(-bore_mm - 6, bore_mm + 6)
    ax.set_ylim(-bore_mm - 8, bore_mm + 8)
    ax.axis("off")
    fig.suptitle(
        "개구부 립 상세 (현재 형상) — 흡입측 OD 평면부 34.8°→13.4° + "
        "R1.5×102.2° + 직선 + R1.0 홈 블렌드 (θ = 0)",
        fontsize=12,
    )

    ax.fill(
        [p[0] / MM for p in contour.material],
        [p[1] / MM for p in contour.material],
        facecolor=ROTOR_COLOR,
        edgecolor="none",
        zorder=2,
    )
    for side in (1.0, -1.0):
        bush = SwingBush()
        xs, ys = _bush_outline_mm(
            state.cutout_center_m[0] + side * bush.piece_shift_m, state.cutout_center_m[1], side
        )
        ax.fill(xs, ys, facecolor=BUSH_COLOR, edgecolor="black", linewidth=0.8, zorder=3)
    vane_x, vane_y = _vane_outline_mm(geometry, state.vane_tip_m[1])
    ax.fill(vane_x, vane_y, facecolor="white", edgecolor="black", linewidth=1.3, zorder=4)
    ax.add_patch(
        plt.Circle(
            (0.0, 0.0), bore_mm, facecolor="none", edgecolor="black", linewidth=2.0, zorder=5
        )
    )
    for edge in (contour.od_arc, contour.inlet_flat, contour.mouth_path):
        ax.plot(
            [p[0] / MM for p in edge],
            [p[1] / MM for p in edge],
            color="black",
            linewidth=2.0,
            zorder=5,
        )

    # Dotted rays from the groove centre (H) to each groove-blend touch point.
    ax.plot(
        groove_x, groove_y, marker="+", color="#c0392b", markersize=10, markeredgewidth=2, zorder=7
    )
    ax.text(
        groove_x + 1.2,
        groove_y,
        "H",
        color="#c0392b",
        fontsize=11,
        ha="left",
        va="center",
        zorder=7,
    )
    for touch in (inlet_touch, outlet_touch):
        ax.plot(
            [groove_x, touch[0]],
            [groove_y, touch[1]],
            color="#2e7d32",
            linestyle=(0, (2, 2)),
            linewidth=1.0,
            zorder=6,
        )
    ax.annotate(
        f"홈 블렌드: 홈 중심 기준 {abs(inlet_angle):.1f}°\n(양측 대칭에 가까움)",
        xy=(inlet_touch[0], inlet_touch[1]),
        xytext=(bore_mm * 0.30, -bore_mm * 0.45),
        fontsize=9,
        color="#2e7d32",
        arrowprops={"arrowstyle": "->", "color": "#2e7d32", "lw": 1.0},
    )

    for (cx_m, cy_m), text in (((0.0, 0.0), "C"), (state.rotor_center_m, "R")):
        cx, cy = cx_m / MM, cy_m / MM
        ax.plot(cx, cy, marker="+", color=GUIDE_COLOR, markersize=9, markeredgewidth=2, zorder=7)
        ax.text(
            cx + 1.4, cy, text, color=GUIDE_COLOR, fontsize=10, ha="left", va="center", zorder=7
        )

    # Zoom box around the inlet lip.
    lip = mouth[: _INLET_LIP_LEN + 3]
    pad = 1.2
    x0, x1 = min(p[0] for p in lip) - pad, max(p[0] for p in lip) + pad
    y0, y1 = min(p[1] for p in lip) - pad, max(p[1] for p in lip) + pad
    ax.add_patch(
        Rectangle(
            (x0, y0),
            x1 - x0,
            y1 - y0,
            facecolor="none",
            edgecolor="#888",
            linewidth=1.0,
            linestyle=(0, (4, 3)),
            zorder=8,
        )
    )

    # Inset: the inlet-lip G1 chain, coloured by segment.
    inset = fig.add_axes([0.62, 0.06, 0.36, 0.52])
    inset.set_aspect("equal")
    inset.set_xlim(x0 - 3.8, x1 + 0.4)
    inset.set_ylim(y0, y1)
    inset.set_xticks([])
    inset.set_yticks([])
    inset.set_title("흡입측 립 확대 — R1.5×102.2° + 직선 + R1.0 (G1 접선 연속)", fontsize=9)
    for edge in (contour.od_arc, contour.inlet_flat):
        inset.plot(
            [p[0] / MM for p in edge], [p[1] / MM for p in edge], color="black", linewidth=1.8
        )
    arc = mouth[1:42]
    straight = mouth[41:43]
    blend = mouth[42:_INLET_LIP_LEN]
    inset.plot([p[0] for p in arc], [p[1] for p in arc], color="#c0392b", linewidth=2.4)
    inset.plot([p[0] for p in straight], [p[1] for p in straight], color="black", linewidth=2.4)
    inset.plot([p[0] for p in blend], [p[1] for p in blend], color="#2e7d32", linewidth=2.4)
    inset.plot(*mouth[0], "o", color="black", markersize=5)
    inset.plot(*mouth[42], "o", color="#2e7d32", markersize=5)
    inset.plot(*inlet_touch, "o", color="#2e7d32", markersize=5)
    inset.annotate(
        "OD 평면부 끝\n(13.4°)",
        xy=mouth[0],
        xytext=(mouth[0][0] - 2.4, mouth[0][1] + 0.4),
        fontsize=8,
        color="black",
    )
    inset.text(
        arc[15][0] - 2.6,
        arc[15][1],
        "R1.5\n(102.2°)",
        fontsize=8,
        color="#c0392b",
        ha="right",
        va="center",
    )
    inset.text(
        straight[0][0] - 2.2,
        straight[0][1] - 0.3,
        "직선 0.40",
        fontsize=8,
        color="black",
        ha="right",
    )
    inset.text(
        blend[15][0] + 0.4,
        blend[15][1],
        "R1.0\n(52.6°)",
        fontsize=8,
        color="#2e7d32",
        ha="left",
        va="center",
    )
    inset.annotate(
        f"D (홈, {abs(inlet_angle):.1f}°)",
        xy=inlet_touch,
        xytext=(inlet_touch[0] + 0.6, inlet_touch[1] - 1.6),
        fontsize=8,
        color="#2e7d32",
        arrowprops={"arrowstyle": "->", "color": "#2e7d32", "lw": 0.9},
    )

    fig.savefig(path)
    plt.close(fig)


# --------------------------------------------------------------------------
# Ported legacy figures (schematics and dimensioned drawings)
# --------------------------------------------------------------------------


def render_dimensioned_side_section(geometry: RotaryGeometry, path: Path) -> None:
    """Axial side section at the vane centre plane (x = 0, theta = 0).

    Hatched end plates and cylinder walls, the gray rotor body with its 14.7
    recess channel, the green 4 MPa(abs) sealed bands, the stepped vane
    (full 21 near the bore, then 2.4 ledges), a red-dashed detail callout and
    an inset detailing the 2.4 / 3.15 / 0.75 stack. All lengths in mm.
    """

    # ---- dimensions straight from the model (mm) -----------------------
    bore_r = geometry.cylinder_radius_m / MM  # 38.5
    ecc = geometry.eccentricity_m / MM  # 4.5
    H = geometry.cylinder_height_m / MM  # 21.0
    hh = 0.5 * H  # 10.5
    bore_d = 2.0 * bore_r  # 77.0

    bands = AxialBands()
    full_depth = bands.full_vane_depth_m / MM  # 15.4
    ledge = bands.ledge_thickness_m / MM  # 2.4
    gap16 = H - 2.0 * ledge  # 16.2
    recess = 3.15  # per side (PHYSICS 3.3)
    tab14 = H - 2.0 * recess  # 14.7
    clear = recess - ledge  # 0.75
    tip_from_bore = 25.0

    # ---- x layout: left inner wall = 0, right inner wall (bore) = 77 ---
    x_wallL = 0.0
    x_rotorL = ecc * 2.0  # 9.0 crescent
    x_R = bore_r + ecc  # 43.0 rotor centre
    x_tip = bore_d - tip_from_bore  # 52.0 vane tip
    x_fulldepth = bore_d - full_depth  # 61.6
    x_tabR = 59.5  # recessed-tab right edge
    x_boreR = bore_d  # 77.0
    hh_tab = 0.5 * tab14  # 7.35
    hh_gap = 0.5 * gap16  # 8.1

    # ---- colours -------------------------------------------------------
    HATCH_FACE = "#eef0f3"
    TAB_GRAY = "#c0c0c0"
    GREEN = "#86c28a"
    DIM = "#20486e"
    RED = "#c0392b"
    GUIDE = "#b7bcc4"
    plate_t = 6.0
    wall_t = 7.0
    x_plateL = x_wallL - wall_t
    x_outerR = x_boreR + 1.4
    x_plateR = x_outerR + 3.0

    plt.rcParams["hatch.linewidth"] = 1.0
    fig = plt.figure(figsize=(13.0, 8.4), dpi=180)
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_aspect("equal")
    ax.set_xlim(-15.0, 92.0)
    ax.set_ylim(-32.0, 33.0)
    ax.axis("off")

    def hatched(x, y, w, h, z=3):
        ax.add_patch(
            Rectangle(
                (x, y),
                w,
                h,
                facecolor=HATCH_FACE,
                edgecolor="black",
                hatch="////",
                linewidth=1.3,
                zorder=z,
            )
        )

    def solid(x, y, w, h, color, z=2, edge="black", lw=1.3):
        ax.add_patch(
            Rectangle((x, y), w, h, facecolor=color, edgecolor=edge, linewidth=lw, zorder=z)
        )

    # full-height faint guide lines (vane tip / full-depth / bore)
    for gx, y0, y1 in ((x_tip, -21.0, 26.0), (x_fulldepth, -12.0, 22.5), (x_boreR, -31.0, 26.0)):
        ax.plot([gx, gx], [y0, y1], color=GUIDE, lw=0.8, zorder=0)

    # ---- green sealed spaces (draw first) ------------------------------
    solid(x_R, hh_tab, x_tip - x_R, hh - hh_tab, GREEN, z=1, edge="none")  # top recess
    solid(x_R, -hh, x_tip - x_R, hh - hh_tab, GREEN, z=1, edge="none")  # bottom recess
    solid(x_tip, hh_tab, x_tabR - x_tip, clear, GREEN, z=1, edge="none")  # top clearance
    solid(x_tip, -hh_gap, x_tabR - x_tip, clear, GREEN, z=1, edge="none")  # bottom clearance

    # ---- rotor body ----------------------------------------------------
    solid(x_rotorL, -hh, x_R - x_rotorL, H, ROTOR_COLOR, z=2)  # main block (21)
    ax.text(
        0.5 * (x_rotorL + x_R),
        0.0,
        "로터 (21 부분)",
        color="#333333",
        ha="center",
        va="center",
        fontsize=14,
        zorder=5,
    )
    solid(x_R, -hh_tab, x_tabR - x_R, tab14, TAB_GRAY, z=2)  # recessed tab (14.7)

    # ---- hatched fixed structure --------------------------------------
    hatched(x_plateL, hh, x_plateR - x_plateL, plate_t)  # top plate
    hatched(x_plateL, -hh - plate_t, x_plateR - x_plateL, plate_t)  # bottom plate
    hatched(x_plateL, -hh, wall_t, H)  # left wall
    hatched(x_outerR, -hh, x_plateR - x_outerR, H)  # outer right wall
    hatched(x_fulldepth, -hh, x_boreR - x_fulldepth, H)  # vane full block
    hatched(x_tip, hh - ledge, x_fulldepth - x_tip, ledge)  # top ledge
    hatched(x_tip, -hh, x_fulldepth - x_tip, ledge)  # bottom ledge

    # dashed section edge on the recessed tab
    ax.plot(
        [x_tabR, x_tabR], [-hh_tab, hh_tab], color="black", lw=1.0, linestyle=(0, (3, 2)), zorder=4
    )

    # ---- dimension helpers --------------------------------------------
    def arrow(x0, y0, x1, y1, color=DIM):
        ax.annotate(
            "",
            xy=(x1, y1),
            xytext=(x0, y0),
            zorder=6,
            arrowprops=dict(arrowstyle="<->", color=color, lw=1.3, shrinkA=0, shrinkB=0),
        )

    def ext(x0, y0, x1, y1):
        ax.plot([x0, x1], [y0, y1], color=GUIDE, lw=0.8, zorder=0)

    # 21 (height, left)
    arrow(-10.0, -hh, -10.0, hh)
    ax.text(-11.6, 0.0, "21", color=DIM, rotation=90, ha="center", va="center", fontsize=13)

    # 25.0 and 15.4 (top right)
    ext(x_tip, hh + plate_t, x_tip, 25.6)
    ext(x_boreR, hh + plate_t, x_boreR, 25.6)
    ext(x_fulldepth, hh + plate_t, x_fulldepth, 22.1)
    arrow(x_tip, 25.0, x_boreR, 25.0)
    ax.text(0.5 * (x_tip + x_boreR), 26.3, "25.0", color=DIM, ha="center", va="bottom", fontsize=13)
    arrow(x_fulldepth, 21.5, x_boreR, 21.5)
    ax.text(
        0.5 * (x_fulldepth + x_boreR),
        22.3,
        "15.4",
        color=DIM,
        ha="center",
        va="bottom",
        fontsize=13,
    )

    # 14.7 and 16.2 (inner, vertical)
    arrow(49.0, -hh_tab, 49.0, hh_tab)
    ax.text(47.7, 0.0, "14.7", color=DIM, rotation=90, ha="center", va="center", fontsize=12)
    arrow(57.0, -hh_gap, 57.0, hh_gap)
    ax.text(55.7, 0.0, "16.2", color=DIM, rotation=90, ha="center", va="center", fontsize=12)

    # a_v = 9 (red, bottom)
    ext(x_R, -hh - plate_t, x_R, -18.6)
    ext(x_tip, -18.6, x_tip, -17.3)
    arrow(x_R, -19.0, x_tip, -19.0, color=RED)
    ax.text(x_R, -17.4, "R", color=RED, ha="center", va="bottom", fontsize=13, fontweight="bold")
    ax.text(49.5, -20.7, "a_v = 9 (θ=0, R→팁)", color=RED, ha="center", va="center", fontsize=12)

    # 68 (rotor OD width)
    ext(x_rotorL, -hh - plate_t, x_rotorL, -23.9)
    ext(x_boreR, -19.4, x_boreR, -23.9)
    arrow(x_rotorL, -23.4, x_boreR, -23.4)
    ax.text(34.0, -22.9, "68 (로터 OD 폭)", color=DIM, ha="center", va="bottom", fontsize=13)

    # 77 (cylinder bore)
    ext(x_wallL, -hh - plate_t, x_wallL, -27.7)
    ext(x_boreR, -23.9, x_boreR, -27.7)
    arrow(x_wallL, -27.2, x_boreR, -27.2)
    ax.text(33.0, -29.9, "77 (실린더 내경)", color=DIM, ha="center", va="top", fontsize=13)

    # ---- detail callout circle + label --------------------------------
    cx, cy, cr = 59.0, 7.2, 6.0
    ax.add_patch(
        plt.Circle(
            (cx, cy), cr, facecolor="none", edgecolor=RED, lw=1.6, linestyle=(0, (5, 4)), zorder=6
        )
    )
    ax.annotate(
        "상세 A",
        xy=(cx - 2.0, cy + cr - 0.5),
        xytext=(50.0, 20.5),
        color=RED,
        fontsize=14,
        ha="center",
        va="bottom",
        arrowprops=dict(arrowstyle="-", color=RED, lw=1.2),
        zorder=6,
    )

    # ---- legend swatch -------------------------------------------------
    solid(-13.0, -20.0, 5.0, 2.6, GREEN, z=2, edge="none")
    ax.text(
        -6.6, -18.7, "4 MPa(절대) 밀폐 공간", color="#333333", ha="left", va="center", fontsize=13
    )

    # top-right caption + title
    ax.text(82.0, 16.0, "상·하 단판", color="#333333", ha="left", va="center", fontsize=13)
    ax.text(
        0.5 * (x_plateL + x_boreR),
        30.5,
        "측면 단면도 — 베인 중심면 x = 0, θ = 0 (단위: mm)",
        color="#111111",
        ha="center",
        va="center",
        fontsize=16,
    )

    # ================= detail A inset ==================================
    iax = fig.add_axes([0.705, 0.02, 0.285, 0.245])
    iax.set_aspect("equal")
    iax.set_xlim(43.0, 78.0)
    iax.set_ylim(2.5, 17.5)
    iax.set_xticks([])
    iax.set_yticks([])
    for s in iax.spines.values():
        s.set_edgecolor("black")
        s.set_linewidth(1.2)
    iax.set_title("상세 A — 턱 2.4 / 리세스 3.15 / 유격 0.75", color=RED, fontsize=13)

    def ihatch(x, y, w, h):
        iax.add_patch(
            Rectangle(
                (x, y),
                w,
                h,
                facecolor=HATCH_FACE,
                edgecolor="black",
                hatch="////",
                linewidth=1.3,
                zorder=3,
            )
        )

    # green recess + clearance
    iax.add_patch(
        Rectangle(
            (x_R, hh_tab), x_tip - x_R, hh - hh_tab, facecolor=GREEN, edgecolor="none", zorder=1
        )
    )
    iax.add_patch(
        Rectangle(
            (x_tip, hh_tab), x_tabR - x_tip, clear, facecolor=GREEN, edgecolor="none", zorder=1
        )
    )
    # rotor tab (gray) below
    iax.add_patch(
        Rectangle(
            (x_R, hh_tab - 3.2),
            x_tabR - x_R,
            3.2,
            facecolor=TAB_GRAY,
            edgecolor="black",
            linewidth=1.2,
            zorder=2,
        )
    )
    # fixed structure
    ihatch(x_plateL, hh, x_plateR - x_plateL, plate_t)  # top plate
    ihatch(x_tip, hh - ledge, x_fulldepth - x_tip, ledge)  # top ledge
    ihatch(x_fulldepth, hh_tab - 3.2, x_boreR - x_fulldepth, hh - (hh_tab - 3.2))  # vane full

    # inset dimensions
    def iarrow(x0, y0, x1, y1, color=DIM):
        iax.annotate(
            "",
            xy=(x1, y1),
            xytext=(x0, y0),
            zorder=6,
            arrowprops=dict(arrowstyle="<->", color=color, lw=1.3, shrinkA=0, shrinkB=0),
        )

    iarrow(46.0, hh_tab, 46.0, hh)
    iax.text(
        44.8,
        0.5 * (hh_tab + hh),
        "3.15",
        color=DIM,
        rotation=90,
        ha="center",
        va="center",
        fontsize=12,
    )
    iarrow(x_fulldepth + 1.2, hh - ledge, x_fulldepth + 1.2, hh)
    iax.text(
        x_fulldepth + 2.2,
        hh - 0.5 * ledge,
        "2.4",
        color=DIM,
        rotation=90,
        ha="center",
        va="center",
        fontsize=12,
    )
    iarrow(55.5, hh_tab, 55.5, hh_tab + clear)
    iax.text(56.5, hh_tab + 0.5 * clear, "0.75", color=DIM, ha="left", va="center", fontsize=12)

    fig.savefig(path)
    plt.close(fig)


def render_stepped_vane_structure(geometry: RotaryGeometry, path: Path) -> None:
    """Stepped vane + thin-rotor channel: axial side section and top view (theta = 0).

    The left panel is a schematic axial section through the vane centre plane
    (x = 0): hatched end plates and cylinder, the full-thickness rotor and its
    14.7 mm thin channel, the 3.15 mm 4 MPa recess pockets (green), and the
    2.4 mm vane ledges. The right panel is the true top-view geometry -- rotor
    contour, both swing-bush pieces, the dark stepped vane, the dotted
    ledge-only zone, and the 8*sqrt(2) = 11.31 mm recess channel.
    """

    _use_korean_font()

    seal_color = "#86c28a"  # 4 MPa sealed recess pockets
    thin_color = "#c0c0c0"  # thin rotor core (박육부)
    plate_face = "#eef0f3"  # hatched-plate fill
    dim_color = "#2c3e50"  # dimension arrows
    red = "#cc0000"

    bands = AxialBands()
    bore_mm = geometry.cylinder_radius_m / MM
    half_h = 0.5 * geometry.cylinder_height_m / MM
    full_depth = bands.full_vane_depth_m / MM
    ledge = bands.ledge_thickness_m / MM
    core_half = 0.5 * 14.7  # thin-rotor core half-height (14.7 mm channel)
    recess_h = half_h - core_half  # 3.15 mm 4 MPa recess per side

    fig = plt.figure(figsize=(16.4, 8.7), dpi=145)
    fig.suptitle(
        "이해 확인 — 단차(스텝) 베인 + 로터 박육 채널 구조: 간섭 없음, "
        "a_v = 9 mm는 턱 팁 기준 (모든 부위 실치수 비례)",
        fontsize=15,
    )

    # ----------------------------------------------------------- left: section
    axL = fig.add_axes([0.015, 0.06, 0.47, 0.80])
    axL.set_aspect("equal")
    axL.axis("off")
    axL.set_xlim(-49, 57)
    axL.set_ylim(-27, 31)
    fig.text(
        0.245,
        0.855,
        "측면 단면 (베인 중심면 x = 0, θ = 0) — 단차 베인 + 박육 로터",
        ha="center",
        va="center",
        fontsize=12.5,
    )

    Y_bore = bore_mm
    Y_full = bore_mm - full_depth
    Y_tip = bore_mm - 25.0
    Y_R = geometry.eccentricity_m / MM
    Y_rbot = Y_R - geometry.rotor_radius_m / MM
    Y_bwall = -bore_mm
    plate_t = 6.0

    def hrect(x, y, w, h):
        axL.add_patch(
            Rectangle(
                (x, y), w, h, facecolor=plate_face, edgecolor="black", linewidth=1.1, hatch="////"
            )
        )

    hrect(-47, half_h, 101, plate_t)
    hrect(-47, -half_h - plate_t, 101, plate_t)
    hrect(Y_bwall - 6, -half_h, 6, 2 * half_h)
    hrect(Y_full, -half_h, 54 - Y_full, 2 * half_h)
    axL.plot([Y_bore, Y_bore], [-half_h, half_h], color="black", lw=1.0)

    axL.add_patch(
        Rectangle(
            (Y_rbot, -half_h),
            Y_R - Y_rbot,
            2 * half_h,
            facecolor=ROTOR_COLOR,
            edgecolor="black",
            linewidth=1.4,
        )
    )
    axL.add_patch(
        Rectangle(
            (Y_R, -core_half),
            Y_full - Y_R,
            2 * core_half,
            facecolor=thin_color,
            edgecolor="black",
            linewidth=1.2,
        )
    )
    for zc in (core_half, -half_h):
        axL.add_patch(
            Rectangle(
                (Y_R, zc),
                Y_tip - Y_R,
                recess_h,
                facecolor=seal_color,
                edgecolor="#3f7a43",
                linewidth=1.0,
            )
        )
    for zc in (half_h - ledge, -half_h):
        axL.add_patch(
            Rectangle(
                (Y_tip, zc),
                Y_full - Y_tip,
                ledge,
                facecolor=plate_face,
                edgecolor="black",
                linewidth=1.0,
                hatch="////",
            )
        )

    axL.text(
        Y_R - 16,
        0,
        "로터 21 mm 부분\n(전체 두께)",
        ha="center",
        va="center",
        fontsize=11,
        color="#555",
    )
    axL.text(
        (Y_R + Y_full) / 2, 0, "박육부\n14.7", ha="center", va="center", fontsize=10.5, color="#444"
    )

    axL.annotate(
        "",
        xy=(Y_bore, 22.5),
        xytext=(Y_tip, 22.5),
        arrowprops=dict(arrowstyle="<->", color=dim_color, lw=1.4),
    )
    axL.text(
        (Y_bore + Y_tip) / 2, 23.2, "25.0", ha="center", va="bottom", fontsize=11, color=dim_color
    )
    axL.annotate(
        "",
        xy=(Y_bore, 19.0),
        xytext=(Y_full, 19.0),
        arrowprops=dict(arrowstyle="<->", color=dim_color, lw=1.4),
    )
    axL.text(
        (Y_bore + Y_full) / 2, 19.6, "15.4", ha="center", va="bottom", fontsize=11, color=dim_color
    )

    axL.annotate(
        "4 MPa 고정 공간 (상·하)\n로터 21 mm부·스윙 부시·베인이 밀폐\n(유막 누출은 추후 고려)",
        xy=(Y_R + 4, half_h - recess_h + 0.2),
        xytext=(-48, 26.5),
        fontsize=10.5,
        color="#2f7a34",
        ha="left",
        arrowprops=dict(arrowstyle="->", color="#2f7a34", lw=1.2),
    )
    axL.annotate(
        "베인 전체 두께(21 mm) 구간:\n"
        "보어에서 15.4 mm까지 — 항상 원형 홈\n"
        "구간 안에만 있어 로터와 간섭 없음",
        xy=(Y_full + 3, half_h - 0.3),
        xytext=(1, 27.0),
        fontsize=10.5,
        color="black",
        ha="left",
        arrowprops=dict(arrowstyle="-", color="black", lw=1.0),
    )
    axL.annotate(
        "상·하 2.4 mm 턱 구간(15.4~25 mm):\n"
        "로터 박육부(14.7)의 리세스(3.15/측)\n"
        "속을 지남 → 편측 유격 0.75 mm",
        xy=(Y_R + 6, -half_h + recess_h + 0.2),
        xytext=(-48, -24.0),
        fontsize=10.5,
        color="black",
        ha="left",
        arrowprops=dict(arrowstyle="->", color="black", lw=1.0),
    )
    axL.annotate(
        "스윙 부시(원형 홈 Ø16 안,\n세부 치수 미정)",
        xy=(Y_bore + 4, -half_h + 1.0),
        xytext=(45, -21),
        fontsize=10.5,
        color="#666",
        ha="left",
        arrowprops=dict(arrowstyle="->", color="#666", lw=1.0),
    )

    for xx in (Y_R, Y_tip):
        axL.plot([xx, xx], [-core_half, -18.0], color=red, lw=1.0, ls=(0, (2, 2)))
    axL.annotate(
        "",
        xy=(Y_tip, -18.0),
        xytext=(Y_R, -18.0),
        arrowprops=dict(arrowstyle="<->", color=red, lw=1.4),
    )
    axL.text(
        Y_R - 0.6, -17.2, "R", ha="center", va="bottom", fontsize=11, color=red, fontweight="bold"
    )
    axL.text(
        (Y_R + Y_tip) / 2, -20.0, "a_v = 9 mm", ha="center", va="center", fontsize=11, color=red
    )
    axL.text(
        (Y_R + Y_tip) / 2,
        -21.8,
        "(로터중심 R → 턱 팁, θ=0)",
        ha="center",
        va="center",
        fontsize=10,
        color=red,
    )

    # -------------------------------------------------------------- right: top
    axR = fig.add_axes([0.54, 0.05, 0.44, 0.80])
    axR.set_aspect("equal")
    axR.axis("off")
    axR.set_xlim(-bore_mm - 8, bore_mm + 8)
    axR.set_ylim(-bore_mm - 8, bore_mm + 10)
    axR.set_title("평면(윗면) 뷰 (θ = 0) — 로터 박육 채널과 단차 베인", fontsize=12.5, pad=10)

    state = prescribed_state(geometry, 0.0)
    contour = rotor_contour(geometry, 0.0)
    groove_x, groove_y = (c / MM for c in state.cutout_center_m)
    rotor_cy = state.rotor_center_m[1] / MM
    vane_tip_y = state.vane_tip_m[1] / MM
    full_depth_y = bore_mm - full_depth
    chan_half = 0.5 * geometry.vane_width_m * sqrt(2.0) / MM
    vane_half = 0.5 * geometry.vane_width_m / MM

    axR.fill(
        [p[0] / MM for p in contour.material],
        [p[1] / MM for p in contour.material],
        facecolor=ROTOR_COLOR,
        edgecolor="none",
        zorder=1,
    )
    axR.add_patch(
        Rectangle(
            (-chan_half, rotor_cy),
            2 * chan_half,
            groove_y - rotor_cy,
            facecolor="#cfcfcf",
            edgecolor="#333",
            linewidth=1.2,
            linestyle=(0, (4, 3)),
            zorder=2,
        )
    )
    bush = SwingBush()
    for side in (1.0, -1.0):
        xs, ys = _bush_outline_mm(
            state.cutout_center_m[0] + side * bush.piece_shift_m, state.cutout_center_m[1], side
        )
        axR.fill(xs, ys, facecolor=BUSH_COLOR, edgecolor="black", linewidth=1.0, zorder=3)
    vx, vy = _vane_outline_mm(geometry, state.vane_tip_m[1])
    axR.fill(vx, vy, facecolor=VANE_COLOR, edgecolor="black", linewidth=1.2, zorder=4)
    axR.add_patch(
        Rectangle(
            (-vane_half, vane_tip_y),
            2 * vane_half,
            full_depth_y - vane_tip_y,
            facecolor="#b7b7b7",
            edgecolor="black",
            linewidth=1.0,
            hatch="....",
            zorder=5,
        )
    )
    axR.add_patch(
        plt.Circle((0, 0), bore_mm, facecolor="none", edgecolor="black", linewidth=2.4, zorder=6)
    )
    for cx, cy in ((0.0, rotor_cy), (0.0, 0.0)):
        axR.plot(cx, cy, marker="+", color="black", markersize=11, markeredgewidth=2, zorder=7)

    for ang, txt in ((-27.0, "OUT"), (27.0, "IN")):
        t0 = port_position(geometry.cylinder_radius_m - 1.5 * MM, ang)
        t1 = port_position(geometry.cylinder_radius_m + 2.5 * MM, ang)
        axR.plot(
            [t0[0] / MM, t1[0] / MM], [t0[1] / MM, t1[1] / MM], color="black", lw=2.2, zorder=7
        )
        lp = port_position(geometry.cylinder_radius_m + 6.0 * MM, ang)
        axR.text(
            lp[0] / MM, lp[1] / MM, txt, ha="center", va="center", fontsize=12, fontweight="bold"
        )

    axR.annotate(
        "베인 전체두께 구간\n(진회색, 보어~15.4)",
        xy=(-vane_half, bore_mm - 6),
        xytext=(-bore_mm - 4, bore_mm + 4),
        fontsize=10.5,
        ha="left",
        arrowprops=dict(arrowstyle="->", color="black", lw=1.0),
    )
    axR.annotate(
        "턱 전용 구간(점 무늬,\n15.4~25, 상·하 2.4 mm만)",
        xy=(vane_half, (vane_tip_y + full_depth_y) / 2),
        xytext=(bore_mm - 12, bore_mm + 4),
        fontsize=10.5,
        ha="left",
        arrowprops=dict(arrowstyle="->", color="black", lw=1.0),
    )
    axR.annotate(
        "박육 14.7 mm 채널\n원형 홈 중심 → 로터 중심,\n폭 8√2 = 약  11.31 mm",
        xy=(-chan_half + 1, rotor_cy + 8),
        xytext=(-bore_mm - 6, -bore_mm * 0.45),
        fontsize=10.5,
        ha="left",
        arrowprops=dict(arrowstyle="->", color="black", lw=1.0),
    )
    axR.annotate(
        "스윙 부시(원형 홈 안)",
        xy=(chan_half + 1, groove_y - 3),
        xytext=(bore_mm * 0.30, -bore_mm * 0.7),
        fontsize=10.5,
        color="#777",
        ha="left",
        arrowprops=dict(arrowstyle="->", color="#777", lw=1.0),
    )

    fig.savefig(path)
    plt.close(fig)


def render_dimensioned_top_view(geometry: RotaryGeometry, path: Path) -> None:
    """Dimensioned top view at crank angle 0 (mm), current asymmetric geometry.

    Reconstructs the legacy ``dimensioned_top_view`` figure straight from the
    live model: bore + rotor OD circles, the ``rotor_contour`` mouth, both swing
    bush pieces, the dark vane, the thin-wall 14.7 channel, and the full set of
    dimensions, port ticks, and the current lip / swing-bush spec block.
    """

    NAVY = "#2e4d7b"
    RED = "#c0392b"
    LEADER = "#7a7f88"

    bore_mm = geometry.cylinder_radius_m / MM
    vane_half = 0.5 * geometry.vane_width_m / MM

    state = prescribed_state(geometry, 0.0)
    contour = rotor_contour(geometry, 0.0)
    rotor_cx = state.rotor_center_m[0] / MM
    rotor_cy = state.rotor_center_m[1] / MM

    channel_w = 8.0 * sqrt(2.0)  # 11.31 mm channel width
    outer_half = 0.5 * 14.7  # thin-wall 14.7 slot envelope

    fig, ax = plt.subplots(figsize=(12.0, 11.8), dpi=120)
    ax.set_aspect("equal")
    ax.set_xlim(-52.0, 62.0)
    ax.set_ylim(-46.0, 64.0)
    ax.axis("off")
    fig.suptitle("상면도 — θ = 0 (단위: mm)", fontsize=15, y=0.965)

    # Rotor material fill (gray) under the channel and solids.
    ax.fill(
        [p[0] / MM for p in contour.material],
        [p[1] / MM for p in contour.material],
        facecolor=ROTOR_COLOR,
        edgecolor="none",
        zorder=1.5,
    )

    # Thin-wall slot envelope (dashed) + the 14.7 dotted channel. The envelope
    # top is tucked up behind the bush so its top edge does not cross the vane.
    ax.add_patch(
        Rectangle(
            (-outer_half, -9.0),
            2.0 * outer_half,
            35.0,
            facecolor="#ececec",
            edgecolor="#9aa0a8",
            linewidth=1.1,
            linestyle=(0, (5, 4)),
            zorder=2.0,
        )
    )
    ax.add_patch(
        Rectangle(
            (-0.5 * channel_w, -2.0),
            channel_w,
            14.7,
            facecolor="#d9d9d9",
            edgecolor="#5c6068",
            hatch="..",
            linewidth=1.0,
            linestyle=(0, (4, 3)),
            zorder=2.2,
        )
    )

    # Swing bush pieces (blue-tinted), then the dark vane over them.
    bush = SwingBush()
    for side in (1.0, -1.0):
        xs, ys = _bush_outline_mm(
            state.cutout_center_m[0] + side * bush.piece_shift_m,
            state.cutout_center_m[1],
            side,
        )
        ax.fill(xs, ys, facecolor=BUSH_COLOR, edgecolor="black", linewidth=0.9, zorder=3)
    vane_x, vane_y = _vane_outline_mm(geometry, state.vane_tip_m[1])
    ax.fill(vane_x, vane_y, facecolor=VANE_COLOR, edgecolor="black", linewidth=1.2, zorder=4)

    # Bore + rotor OD outline (two black circles) and the asymmetric mouth.
    ax.add_patch(
        plt.Circle(
            (0.0, 0.0), bore_mm, facecolor="none", edgecolor="black", linewidth=2.4, zorder=5
        )
    )
    for edge in (contour.od_arc, contour.inlet_flat, contour.mouth_path):
        ax.plot(
            [p[0] / MM for p in edge],
            [p[1] / MM for p in edge],
            color="black",
            linewidth=1.8,
            zorder=5,
        )

    # C (bore centre) and R (rotor centre) markers.
    ax.plot(0.0, 0.0, marker="+", color=GUIDE_COLOR, markersize=11, markeredgewidth=2.2, zorder=7)
    ax.text(
        0.0, -2.6, "C", color=GUIDE_COLOR, fontweight="bold", fontsize=11, ha="center", zorder=7
    )
    ax.plot(
        rotor_cx,
        rotor_cy,
        marker="+",
        color=GUIDE_COLOR,
        markersize=11,
        markeredgewidth=2.2,
        zorder=7,
    )
    ax.text(
        rotor_cx + 1.7,
        rotor_cy - 0.2,
        "R",
        color=GUIDE_COLOR,
        fontweight="bold",
        fontsize=11,
        ha="left",
        va="center",
        zorder=7,
    )

    # e = 4.5 (C -> R).
    ax.annotate(
        "",
        xy=(-13.0, 4.5),
        xytext=(-13.0, 0.0),
        arrowprops=dict(arrowstyle="<->", color=NAVY, lw=1.3),
    )
    ax.text(-14.5, 2.25, "e = 4.5", color=NAVY, fontsize=11, ha="right", va="center")

    # L = 25 (R -> groove centre).
    ax.annotate(
        "",
        xy=(-27.0, 29.5),
        xytext=(-27.0, 4.5),
        arrowprops=dict(arrowstyle="<->", color=NAVY, lw=1.3),
    )
    ax.text(-29.5, 17.0, "L = 25", color=NAVY, fontsize=11, ha="right", va="center")

    # Vane width 8 (top), with witness lines down to the flanks.
    y_top_arrow = bore_mm + 4.0
    ax.annotate(
        "",
        xy=(vane_half, y_top_arrow),
        xytext=(-vane_half, y_top_arrow),
        arrowprops=dict(arrowstyle="<->", color=NAVY, lw=1.3),
    )
    ax.text(
        0.0, y_top_arrow + 1.8, "8 (베인 폭)", color=NAVY, fontsize=11, ha="center", va="bottom"
    )
    for sx in (-vane_half, vane_half):
        ax.plot([sx, sx], [bore_mm - 1.0, y_top_arrow], color=NAVY, lw=0.7)

    # Channel width 8√2 = 11.31 (bottom).
    y_cw = -3.5
    ax.annotate(
        "",
        xy=(0.5 * channel_w, y_cw),
        xytext=(-0.5 * channel_w, y_cw),
        arrowprops=dict(arrowstyle="<->", color=NAVY, lw=1.3),
    )
    ax.text(
        0.0, y_cw - 1.8, "8√2 = 11.31 (채널 폭)", color=NAVY, fontsize=10.5, ha="center", va="top"
    )

    # 15.4 (full vane depth) and 25.0 (tip from bore) on the right.
    y_ref = 25.0
    x15 = 19.0
    ax.annotate(
        "",
        xy=(x15, y_ref - 15.4),
        xytext=(x15, y_ref),
        arrowprops=dict(arrowstyle="<->", color=NAVY, lw=1.3),
    )
    ax.text(x15 + 1.4, y_ref - 7.7, "15.4", color=NAVY, fontsize=11, ha="left", va="center")
    x25 = 27.0
    ax.annotate(
        "",
        xy=(x25, y_ref - 25.0),
        xytext=(x25, y_ref),
        arrowprops=dict(arrowstyle="<->", color=NAVY, lw=1.3),
    )
    ax.text(x25 + 1.4, y_ref - 12.5, "25.0", color=NAVY, fontsize=11, ha="left", va="center")

    # OUT / IN port markers with 30 deg ticks.
    for ang, name, ha in ((-30.0, "OUT", "right"), (30.0, "IN", "left")):
        tin = port_position(geometry.cylinder_radius_m - 1.5 * MM, ang)
        tout = port_position(geometry.cylinder_radius_m + 2.5 * MM, ang)
        ax.plot(
            [tin[0] / MM, tout[0] / MM],
            [tin[1] / MM, tout[1] / MM],
            color="black",
            lw=1.6,
            zorder=6,
        )
        dpos = port_position(geometry.cylinder_radius_m - 4.5 * MM, ang)
        ax.text(
            dpos[0] / MM,
            dpos[1] / MM,
            "30°",
            color="black",
            fontsize=10,
            ha="center",
            va="center",
            zorder=6,
        )
        npos = port_position(geometry.cylinder_radius_m + 6.0 * MM, ang)
        ax.text(
            npos[0] / MM,
            npos[1] / MM,
            name,
            color="black",
            fontweight="bold",
            fontsize=12,
            ha=ha,
            va="center",
            zorder=6,
        )

    # Mouth lip angle callouts: outlet/OUT 13 deg (blue), inlet/IN 13.4 deg (red).
    ax.plot([6.5, 12.5], [19.0, 31.5], color=RED, lw=1.0, linestyle=(0, (2, 2)), zorder=6)
    ax.plot([-6.5, -12.5], [19.0, 31.5], color=RED, lw=1.0, linestyle=(0, (2, 2)), zorder=6)
    ax.text(9.5, 17.3, "13.4°", color=RED, fontsize=10.5, ha="center", va="center", zorder=7)
    ax.text(-9.5, 17.3, "13°", color=NAVY, fontsize=10.5, ha="center", va="center", zorder=7)

    # Leader labels.
    ax.annotate(
        "Ø16 (원형 홈)",
        xy=(-6.0, 30.8),
        xytext=(-42.0, 46.0),
        color=GUIDE_COLOR,
        fontsize=11,
        ha="left",
        va="center",
        arrowprops=dict(arrowstyle="-", color=LEADER, lw=1.0),
    )
    ax.annotate(
        "박육 14.7 채널(점선)",
        xy=(-outer_half, 6.0),
        xytext=(-51.0, 5.0),
        color=GUIDE_COLOR,
        fontsize=11,
        ha="left",
        va="center",
        arrowprops=dict(arrowstyle="-", color=LEADER, lw=1.0),
    )
    ax.annotate(
        "Ø77 (실린더 내경)",
        xy=(-27.2, -27.2),
        xytext=(-51.0, -33.0),
        color=NAVY,
        fontsize=11,
        ha="left",
        va="center",
        arrowprops=dict(arrowstyle="->", color=NAVY, lw=1.0),
    )
    ax.annotate(
        "Ø68 (로터 OD)",
        xy=(-21.9, -21.5),
        xytext=(-51.0, -42.0),
        color=NAVY,
        fontsize=11,
        ha="left",
        va="center",
        arrowprops=dict(arrowstyle="->", color=NAVY, lw=1.0),
    )

    # Current lip + swing-bush spec block (top right).
    spec = (
        "흡입측: OD 평면부 34.8°→13.4°(R33.657, 12.57)\n"
        "  + R1.5×102.2° → 직선 0.399 → R1.0×52.62° → 홈 47.74°\n"
        "토출측: OD 접선 13° + R1.5×92.6° → 직선 0.744\n"
        "  → R1.0×52.43° → 홈 47.97°\n"
        "스윙 부시 2조각: 곡면 Ø15.94 × 106.1°\n"
        "평면부 3.99 (접촉 11.94), 사이 R0.5 블렌드\n"
        "조각 중심 ±0.020 → 유막 0.010/0.010"
    )
    ax.text(
        3.0, 62.5, spec, color=NAVY, fontsize=9.0, ha="left", va="top", linespacing=1.55, zorder=8
    )
    ax.plot([15.0, 8.0], [46.5, 34.0], color=NAVY, lw=0.8, zorder=6)

    fig.savefig(path)
    plt.close(fig)


def render_bush_placement_clearances(geometry: RotaryGeometry, path: Path) -> None:
    """Two panels: top-view swing-bush placement + axial-section 4 MPa seal."""

    bush_blue = "#9fb6d4"
    green = "#2e7d32"
    red = "#c0392b"

    bore_mm = geometry.cylinder_radius_m / MM
    state = prescribed_state(geometry, 0.0)
    groove_x, groove_y = (c / MM for c in state.cutout_center_m)
    bush = SwingBush()

    fig = plt.figure(figsize=(19.5, 10.6), dpi=130)
    ax_l = fig.add_axes([0.02, 0.03, 0.46, 0.86])
    ax_r = fig.add_axes([0.52, 0.03, 0.46, 0.86])

    fig.suptitle(
        "스윙 부시 배치와 클리어런스 — 곡면 Ø15.94×106.1° + R0.5 블렌드 · "
        "평면부 3.99 · 중심 ±0.020 → 유막 0.010/0.010 (단위 mm)",
        fontsize=14,
    )

    # ------------------------------------------------------------------
    # LEFT: top view
    # ------------------------------------------------------------------
    ax_l.set_aspect("equal")
    ax_l.set_xlim(-46, 46)
    ax_l.set_ylim(-47, 45)
    ax_l.axis("off")
    ax_l.set_title("상면 (θ = 0) — 스윙 부시 배치와 포획 조건", fontsize=13, pad=8)

    contour = rotor_contour(geometry, 0.0)
    ax_l.fill(
        [p[0] / MM for p in contour.material],
        [p[1] / MM for p in contour.material],
        facecolor=ROTOR_COLOR,
        edgecolor="none",
        zorder=2,
    )

    # dotted channel rectangle (박육부 채널, 폭 14.7)
    ch_w = 14.7
    ax_l.add_patch(
        Rectangle(
            (-ch_w / 2, groove_y - 10.0),
            ch_w,
            10.5,
            facecolor="none",
            edgecolor="#555",
            linewidth=1.1,
            linestyle=(0, (2, 2)),
            zorder=3,
        )
    )
    ax_l.plot(
        0.0,
        groove_y - 6.0,
        marker="+",
        color=GUIDE_COLOR,
        markersize=11,
        markeredgewidth=2,
        zorder=6,
    )

    # swing-bush pieces (light blue override)
    gx_m, gy_m = state.cutout_center_m
    for side in (1.0, -1.0):
        xs, ys = _bush_outline_mm(gx_m + side * bush.piece_shift_m, gy_m, side)
        ax_l.fill(xs, ys, facecolor=bush_blue, edgecolor="black", linewidth=0.9, zorder=4)

    # dark vane
    vane_x, vane_y = _vane_outline_mm(geometry, state.vane_tip_m[1])
    ax_l.fill(vane_x, vane_y, facecolor=VANE_COLOR, edgecolor="black", linewidth=1.2, zorder=5)

    # bore outline + rotor edges
    ax_l.add_patch(
        plt.Circle(
            (0.0, 0.0), bore_mm, facecolor="none", edgecolor="black", linewidth=2.2, zorder=6
        )
    )
    for edge in (contour.od_arc, contour.inlet_flat, contour.mouth_path):
        ax_l.plot(
            [p[0] / MM for p in edge],
            [p[1] / MM for p in edge],
            color="black",
            linewidth=2.2,
            zorder=6,
        )

    # red H at groove centre
    ax_l.plot(
        groove_x, groove_y, marker="+", color=red, markersize=12, markeredgewidth=2.2, zorder=8
    )
    ax_l.text(
        groove_x + 1.2, groove_y + 0.3, "H", color=red, fontsize=12, fontweight="bold", zorder=8
    )

    # green blend-start rays, ±37.07° from vertical, from H upward
    ray_ang = 37.07
    ray_len = 9.5
    for s in (1.0, -1.0):
        ex = groove_x + s * ray_len * sin(radians(ray_ang))
        ey = groove_y + ray_len * cos(radians(ray_ang))
        ax_l.plot(
            [groove_x, ex],
            [groove_y, ey],
            color=green,
            linestyle=(0, (2, 2)),
            linewidth=1.3,
            zorder=7,
        )

    # red film lines (exaggerated) pointing into the bush
    ax_l.plot([groove_x, -18.0], [groove_y - 2.0, -30.0], color=red, linewidth=1.0, zorder=7)
    ax_l.plot([groove_x + 4.5, 8.0], [groove_y - 3.0, -30.0], color=red, linewidth=1.0, zorder=7)

    def note(text, xy, xytext, color, ha="left"):
        ax_l.annotate(
            text,
            xy=xy,
            xycoords="data",
            xytext=xytext,
            textcoords="axes fraction",
            fontsize=9.5,
            color=color,
            ha=ha,
            va="top",
            arrowprops={"arrowstyle": "->", "color": color, "lw": 1.0},
        )

    note(
        "스윙 부시 2조각:\n조각 중심 = H에서 ±0.020 오프셋(TDC)\n"
        "평면은 베인과 정렬, 로터에 대해 ±10.4° 요동\nx방향 총 유격 0.020 (0.020은 그 중앙)",
        xy=(-6.5, groove_y + 1.0),
        xytext=(-0.01, 0.86),
        color="black",
    )
    note(
        "조각 곡면 끝: 37.07°\n(홈 중심 기준, 블렌드 시작)",
        xy=(groove_x + ray_len * sin(radians(ray_ang)), groove_y + ray_len * cos(radians(ray_ang))),
        xytext=(0.60, 0.90),
        color=green,
    )
    note(
        "포획: 조각 두께 3.980 > 마우스 개구 2.573\n(홈 원 위 베인면~립 블렌드점 거리)\n"
        "→ 빠져나갈 수 없음 (각도 여유는 최대 기울기서 1.07°)",
        xy=(groove_x + 6.0, groove_y + 5.0),
        xytext=(0.58, 0.74),
        color=green,
    )
    note(
        "곡면부-홈 유막 0.010 (운전 위치 최소, 과장 표시)\n"
        "부시 조각 OD R7.970 (Ø15.94), 곡면 106.1°",
        xy=(-18.0, -30.0),
        xytext=(-0.01, 0.30),
        color=red,
    )
    note(
        "평면부-베인 유막 0.010 (과장 표시)\n평면부 3.99, 접촉 길이 11.94",
        xy=(8.0, -30.0),
        xytext=(0.50, 0.24),
        color=red,
    )
    note(
        "부시-베인 상대 슬라이딩:\n행정 9 mm, 최대 0.86 m/s (30 Hz)",
        xy=(-24.0, -24.0),
        xytext=(-0.01, 0.13),
        color="black",
    )

    # ------------------------------------------------------------------
    # RIGHT: axial section
    # ------------------------------------------------------------------
    ax_r.set_aspect("equal")
    ax_r.axis("off")
    ax_r.set_title("축방향 단면 (x = +5 mm, θ = 0) — 부시 높이와 4 MPa 밀폐", fontsize=13, pad=8)

    H = 21.0  # chamber height between plates
    plate_h = 5.0
    x0 = 0.0
    total_w = 66.0
    ax_r.set_xlim(-2, total_w + 2)
    ax_r.set_ylim(-plate_h - 12.0, H + plate_h + 14.0)

    # hatched plates (cylinder heads)
    for y in (H, -plate_h):
        ax_r.add_patch(
            Rectangle(
                (x0, y),
                total_w,
                plate_h,
                facecolor="white",
                edgecolor="black",
                linewidth=1.2,
                hatch="////",
                zorder=2,
            )
        )

    # rotor full-thickness block
    rot_x0, rot_w = 2.0, 30.0
    ax_r.add_patch(
        Rectangle(
            (rot_x0, 0.0),
            rot_w,
            H,
            facecolor=ROTOR_COLOR,
            edgecolor="black",
            linewidth=1.1,
            zorder=3,
        )
    )
    ax_r.text(
        rot_x0 + rot_w / 2,
        H / 2,
        "로터(전체 두께 21)",
        ha="center",
        va="center",
        fontsize=10,
        zorder=5,
    )

    # thin-wall (박육부) block, 14.7 tall, centred -> 3.15 gap top and bottom
    thin_h = 14.7
    thin_gap = (H - thin_h) / 2.0  # 3.15
    thin_x0, thin_w = rot_x0 + rot_w, 16.0
    for gy0, gy1 in ((0.0, thin_gap), (H - thin_gap, H)):
        ax_r.add_patch(
            Rectangle(
                (thin_x0, gy0),
                thin_w + 4.0,
                gy1 - gy0,
                facecolor="#8fce8f",
                edgecolor="none",
                zorder=3,
            )
        )
    ax_r.add_patch(
        Rectangle(
            (thin_x0, thin_gap),
            thin_w,
            thin_h,
            facecolor=ROTOR_COLOR,
            edgecolor="black",
            linewidth=1.1,
            zorder=4,
        )
    )
    ax_r.text(
        thin_x0 + thin_w / 2, H / 2, "박육부 14.7", ha="center", va="center", fontsize=10, zorder=5
    )

    # swing-bush block, height 20.983 -> c_a=0.0085/side (exaggerated gap)
    ca = 0.35
    bush_x0, bush_w = thin_x0 + thin_w + 4.0, 10.0
    ax_r.add_patch(
        Rectangle(
            (bush_x0, ca),
            bush_w,
            H - 2 * ca,
            facecolor=bush_blue,
            edgecolor="black",
            linewidth=1.3,
            zorder=5,
        )
    )
    ax_r.text(
        bush_x0 + bush_w / 2,
        H / 2,
        "스윙 부시 단면",
        ha="center",
        va="center",
        fontsize=10,
        rotation=90,
        zorder=6,
    )

    def rnote(text, xy, xytext, color, ha="left"):
        ax_r.annotate(
            text,
            xy=xy,
            xycoords="data",
            xytext=xytext,
            textcoords="axes fraction",
            fontsize=9.5,
            color=color,
            ha=ha,
            va="top",
            arrowprops={"arrowstyle": "->", "color": color, "lw": 1.0},
        )

    rnote(
        "4 MPa 공간이 부시 OD에서 막힘\n→ 부시를 베인 쪽(+y)으로 미는 배압\n"
        "(홈 둘레 중 채널 개구 ±45° 구간의\n상·하 3.15 mm 띠)",
        xy=(thin_x0 + thin_w, H - thin_gap / 2),
        xytext=(0.30, 0.94),
        color=green,
    )
    rnote(
        "c_a = 0.0085/편측 (과장 표시)\n→ 부시 높이 20.983",
        xy=(bush_x0 + bush_w / 2, H - ca / 2),
        xytext=(0.66, 0.94),
        color=red,
    )
    rnote(
        "밀폐를 위해 부시는 전체 높이(약 21)를\n채워야 함 — 14.7이면 4 MPa 띠가 열림",
        xy=(thin_x0 + thin_w, thin_gap / 2),
        xytext=(0.06, 0.16),
        color="#1a3a6a",
    )
    rnote(
        "베인(전체 두께)-부시 물림 길이:\n"
        "TDC 14.4 mm → BDC 5.4 mm\n"
        "(x=0 단면에서, 그림은 x=+5 단면)",
        xy=(bush_x0 + bush_w, 3.0),
        xytext=(0.62, 0.16),
        color="#333333",
    )

    fig.savefig(path)
    plt.close(fig)


def render_vane_side_view(geometry: RotaryGeometry, path: Path) -> None:
    """Axial-section side view of the vane, vertical-shaft frame.

    Horizontal = y (radial / vane-length direction), vertical = z (axial).
    Panel A is the rejected separate-sliding-vane reading; panel B is the
    confirmed cylinder-integral fixed vane, shown with the contradiction it
    creates when the tip distance a_v is applied.
    """

    RED = "#d92b2b"
    GRAY = GUIDE_COLOR
    ORANGE = COMPRESSION_COLOR

    h_mm = geometry.cylinder_height_m / MM  # 21 axial height
    vane_w = geometry.vane_width_m / MM  # 8 into-page width
    a_v = geometry.vane_tip_distance_at_top_m / MM  # 9 centre-to-tip
    ecc = geometry.eccentricity_m / MM  # 4.5
    tip_y = ecc + a_v  # 13.5
    hole_d = 2.0 * geometry.cutout_radius_m / MM  # 16 (rotor hole diameter)

    # Shared vertical layout: two fixed end plates enclose the axial gap H.
    plate_t = 3.0
    z0, z1 = 0.0, plate_t  # bottom plate
    z2, z3 = plate_t, plate_t + h_mm  # inner gap (H)
    z4 = z3 + plate_t  # top plate top
    x_left, x_right = 0.0, 50.0

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(16.0, 8.2), dpi=120)
    fig.suptitle(
        "베인 측면(축방향 단면) 뷰 — 수직 샤프트 기준, "
        "가로 = y(반경·베인 길이 방향), 세로 = z(축방향)",
        fontsize=14,
    )

    def draw_plates(ax):
        for zb in (z0, z3):
            ax.add_patch(
                Rectangle(
                    (x_left, zb),
                    x_right - x_left,
                    plate_t,
                    facecolor="white",
                    edgecolor="black",
                    hatch="////",
                    linewidth=1.2,
                    zorder=3,
                )
            )

    def draw_h_dim(ax):
        xd = 2.2
        ax.annotate(
            "",
            xy=(xd, z2 + 0.4),
            xytext=(xd, z3 - 0.4),
            arrowprops={"arrowstyle": "<->", "color": "black", "lw": 1.4},
            zorder=6,
        )
        ax.text(
            xd - 1.0,
            (z2 + z3) / 2.0,
            f"H = {h_mm:.0f} mm\n(축방향)",
            ha="right",
            va="center",
            fontsize=10,
        )

    for ax in (ax_a, ax_b):
        ax.set_xlim(-8.5, 57.0)
        ax.set_ylim(-13.5, 45.0)
        ax.set_aspect("equal")
        ax.axis("off")

    # ------------------------------------------------------------------
    # Panel A: separate sliding vane (rejected)
    # ------------------------------------------------------------------
    ax_a.set_title("A. 분리형 슬라이딩 베인 — 측면 단면 (x = 0, θ = 0)", fontsize=12)
    draw_plates(ax_a)
    draw_h_dim(ax_a)

    body_l, body_r = 5.0, 26.0
    band_r = body_r + vane_w * 0.42  # narrow interference band
    vane_l, vane_r = band_r, band_r + 11.0

    ax_a.add_patch(
        Rectangle(
            (body_l, z2),
            body_r - body_l,
            h_mm,
            facecolor=ROTOR_COLOR,
            edgecolor="black",
            linewidth=1.0,
            zorder=2,
        )
    )
    ax_a.text(
        (body_l + body_r) / 2.0,
        (z2 + z3) / 2.0,
        "로터 몸체 단면 (x = 0)",
        ha="center",
        va="center",
        fontsize=10,
    )
    # Red cross-hatched interference band at the tip / rotor-body overlap.
    ax_a.add_patch(
        Rectangle(
            (body_r, z2),
            band_r - body_r,
            h_mm,
            facecolor="white",
            edgecolor=RED,
            hatch="xxxx",
            linewidth=1.4,
            zorder=4,
        )
    )
    # Orange vane plate (full axial height, slides in y).
    ax_a.add_patch(
        Rectangle(
            (vane_l, z2),
            vane_r - vane_l,
            h_mm,
            facecolor=ORANGE,
            edgecolor="black",
            linewidth=1.2,
            zorder=3,
        )
    )
    # Dotted "wall position" lines (through slot -> no wall in this section).
    for xw in (vane_r - 4.0, vane_r - 1.0):
        ax_a.plot([xw, xw], [z2, z3], ls=":", color=GRAY, linewidth=1.2, zorder=5)

    # s(theta) reciprocation arrow.
    ax_a.annotate(
        "",
        xy=(vane_r + 8.5, (z2 + z3) / 2.0),
        xytext=(vane_r + 1.0, (z2 + z3) / 2.0),
        arrowprops={"arrowstyle": "<->", "color": RED, "lw": 1.8},
        zorder=6,
    )
    ax_a.text(
        vane_r + 4.7,
        (z2 + z3) / 2.0 + 3.0,
        "s(θ): y방향\n9 mm 왕복",
        ha="center",
        va="bottom",
        fontsize=9.5,
        color=RED,
    )

    # Notes.
    ax_a.annotate(
        "베인(주황): 두께 8 mm(지면 수직 방향),\n"
        "축방향 전체 높이 21 mm 판 —\n"
        "단판 사이에서 y방향으로만 미끄러짐",
        xy=((vane_l + vane_r) / 2.0, z3),
        xytext=(1.0, z4 + 9.0),
        fontsize=9.5,
        color=ORANGE,
        ha="left",
        va="bottom",
        arrowprops={"arrowstyle": "->", "color": ORANGE, "lw": 1.0},
    )
    ax_a.annotate(
        "벽 위치(점선): 관통 슬롯이라\n이 단면(x=0)에는 벽이 없음",
        xy=(vane_r - 1.0, z3),
        xytext=(27.0, z4 + 11.0),
        fontsize=9.5,
        color=GRAY,
        ha="left",
        va="bottom",
        arrowprops={"arrowstyle": "->", "color": GRAY, "lw": 1.0},
    )
    ax_a.annotate(
        "상·하 단판(빗금): 고정",
        xy=(46.0, z4),
        xytext=(57.0, z4 + 6.0),
        fontsize=9.5,
        color="black",
        ha="right",
        va="bottom",
        arrowprops={"arrowstyle": "->", "color": "black", "lw": 1.0},
    )
    ax_a.annotate(
        f"팁~로터 몸체 간섭 약 {vane_w:.0f} mm (전 각도)\n— 현재 GUI 투명 겹침, CAD 확인 대기",
        xy=((body_r + band_r) / 2.0, z2),
        xytext=(1.0, z1 - 10.0),
        fontsize=9.5,
        color=RED,
        ha="left",
        va="top",
        arrowprops={"arrowstyle": "->", "color": RED, "lw": 1.0},
    )

    # ------------------------------------------------------------------
    # Panel B: cylinder-integral fixed vane (confirmed model)
    # ------------------------------------------------------------------
    ax_b.set_title("B. 실린더 일체 고정 베인 — 측면 단면 (x = 0, θ = 0)", fontsize=12)
    draw_plates(ax_b)
    draw_h_dim(ax_b)

    body_l, body_r = 5.0, 32.0
    hole_edge = 41.0
    rib_l, rib_r = 44.0, 49.5

    ax_b.add_patch(
        Rectangle(
            (body_l, z2),
            body_r - body_l,
            h_mm,
            facecolor=ROTOR_COLOR,
            edgecolor="black",
            linewidth=1.0,
            zorder=2,
        )
    )
    ax_b.text(
        (body_l + body_r) / 2.0,
        (z2 + z3) / 2.0,
        "로터 몸체 단면 (x = 0)",
        ha="center",
        va="center",
        fontsize=10,
    )
    # Integral fixed rib from the bore (hatched, connects both plates).
    ax_b.add_patch(
        Rectangle(
            (rib_l, z2),
            rib_r - rib_l,
            h_mm,
            facecolor="white",
            edgecolor="black",
            hatch="////",
            linewidth=1.2,
            zorder=3,
        )
    )
    # Dotted rotor-hole boundary.
    ax_b.plot([hole_edge, hole_edge], [z2, z3], ls=":", color=GRAY, linewidth=1.3, zorder=5)
    # Red X at the rotor-body edge + dashed line to the hole boundary.
    zc = (z2 + z3) / 2.0
    ax_b.plot([body_r, hole_edge], [zc, zc], ls="--", color=RED, linewidth=1.3, zorder=6)
    ax_b.plot(body_r, zc, marker="x", color=RED, markersize=11, markeredgewidth=3, zorder=7)

    # Rotor-hole reciprocation arrow (below the assembly).
    ax_b.annotate(
        "",
        xy=(46.0, z1 - 6.0),
        xytext=(36.0, z1 - 6.0),
        arrowprops={"arrowstyle": "<->", "color": RED, "lw": 1.8},
        zorder=6,
    )
    ax_b.text(
        41.0,
        z1 - 8.0,
        "θ에 따라 로터 구멍이 y로 9 mm 왕복",
        ha="center",
        va="top",
        fontsize=9.5,
        color=RED,
    )

    # Notes.
    ax_b.annotate(
        "베인: 실린더 보어에서 안쪽으로 뻗은\n"
        f"일체형 리브(빗금) — 고정, 전체 높이 {h_mm:.0f} mm\n"
        "(팁은 항상 구멍 안: 그림 y = 25 mm)",
        xy=(rib_l, z3),
        xytext=(1.0, z4 + 9.0),
        fontsize=9.5,
        color="black",
        ha="left",
        va="bottom",
        arrowprops={"arrowstyle": "->", "color": "black", "lw": 1.0},
    )
    ax_b.annotate(
        f"점선: 로터 구멍(Ø{hole_d:.0f}) 경계",
        xy=(hole_edge, z3),
        xytext=(30.0, z4 + 11.0),
        fontsize=9.5,
        color=GRAY,
        ha="left",
        va="bottom",
        arrowprops={"arrowstyle": "->", "color": GRAY, "lw": 1.0},
    )
    ax_b.annotate(
        "상·하 단판(빗금): 고정",
        xy=(46.0, z4),
        xytext=(57.0, z4 + 6.0),
        fontsize=9.5,
        color="black",
        ha="right",
        va="bottom",
        arrowprops={"arrowstyle": "->", "color": "black", "lw": 1.0},
    )
    ax_b.annotate(
        f"a_v = {a_v:.0f} mm 적용 시 팁(θ=0)은 y = {tip_y:.1f}\n→ 로터 몸체 내부(X), 일체형과 모순",
        xy=(body_r, zc),
        xytext=(1.0, z1 - 10.0),
        fontsize=9.5,
        color=RED,
        ha="left",
        va="top",
        arrowprops={"arrowstyle": "->", "color": RED, "lw": 1.0},
    )

    fig.text(
        0.5,
        0.035,
        "두 해석 모두 베인은 두께 8 mm × 축방향 21 mm 판이며 상·하 단판 사이 전체 높이를 채움.  "
        "A: 판이 단판 사이에서 y방향으로 미끄러짐(별도 부품).  "
        "B: 판이 실린더와 한 몸(고정)이고 로터 구멍이 그 위를 미끄러짐.",
        ha="center",
        va="bottom",
        fontsize=10.5,
    )

    fig.tight_layout(rect=(0.0, 0.07, 1.0, 0.95))
    fig.savefig(path)
    plt.close(fig)


def render_vane_model_comparison(geometry: RotaryGeometry, path: Path) -> None:
    """Two candidate vane/rotor interpretations at theta = 0 (fixed = hatched,
    moving = plain): A separate sliding vane vs. B the confirmed cylinder-integral
    fixed vane. Both are drawn schematically; B carries the y = 13.5 contradiction."""

    RED = "#cc1a1a"
    DARK = "#1a1a1a"
    VANE_EDGE = "#7b3f1e"
    HATCH_BG = "#f2f2f2"

    bore = geometry.cylinder_radius_m / MM  # 38.5
    rotor_r = geometry.rotor_radius_m / MM  # 34.0
    ecc = bore - rotor_r  # 4.5 (rotor pushed up at TDC)
    wall_out = bore + 6.5  # cylinder OD
    hole_r = 8.0  # rotor cutout radius (Ø16)
    hole_cy = 29.0  # cutout centre, mm above bore centre
    half_w = 0.5 * geometry.vane_width_m / MM  # 4.0

    def draw_body(ax):
        # cylinder wall = hatched annulus (fixed); bore centre at origin
        ax.add_patch(
            plt.Circle(
                (0, 0),
                wall_out,
                facecolor=HATCH_BG,
                edgecolor="black",
                linewidth=2.0,
                hatch="///",
                zorder=1,
            )
        )
        ax.add_patch(plt.Circle((0, 0), bore, facecolor="white", edgecolor="none", zorder=2))
        ax.add_patch(
            plt.Circle((0, 0), bore, facecolor="none", edgecolor="black", linewidth=2.0, zorder=3)
        )
        # rotor (moving, plain grey), pushed up by the eccentricity at theta = 0
        ax.add_patch(
            plt.Circle((0, ecc), rotor_r, facecolor=ROTOR_COLOR, edgecolor="none", zorder=4)
        )
        ax.add_patch(
            plt.Circle(
                (0, ecc), rotor_r, facecolor="none", edgecolor="black", linewidth=2.0, zorder=5
            )
        )
        # rotor cutout hole (Ø16)
        ax.add_patch(
            plt.Circle(
                (0, hole_cy), hole_r, facecolor="white", edgecolor="black", linewidth=1.6, zorder=6
            )
        )
        # centre markers: rotor centre (grey +) and bore centre (black +)
        ax.plot(0, ecc, marker="+", color=GUIDE_COLOR, markersize=11, markeredgewidth=2.4, zorder=9)
        ax.plot(0, 0, marker="+", color="black", markersize=11, markeredgewidth=2.4, zorder=9)

    def frame(ax, title):
        ax.set_aspect("equal")
        ax.set_xlim(-88, 82)
        ax.set_ylim(-62, 70)
        ax.axis("off")
        ax.set_title(title, fontsize=13, linespacing=1.35)

    _use_korean_font()
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(17.9, 9.66), dpi=140)

    # ---------------- Panel A: separate sliding vane ----------------
    frame(
        axL,
        "A. 제가 가정했던 구조 — 분리형 슬라이딩 베인\n"
        "(롤링피스톤식처럼 베인이 실린더 슬롯 안에서 미끄러짐)",
    )
    draw_body(axL)
    # orange vane sticking out through the wall slot
    axL.add_patch(
        Rectangle(
            (-half_w, 14.0),
            2 * half_w,
            48.0 - 14.0,
            facecolor=COMPRESSION_COLOR,
            edgecolor=VANE_EDGE,
            linewidth=2.0,
            zorder=8,
        )
    )
    # s(theta) stroke arrow above the vane
    axL.annotate(
        "",
        xy=(0, 57.5),
        xytext=(0, 50.0),
        arrowprops={"arrowstyle": "<->", "color": RED, "lw": 2.2},
        zorder=10,
    )
    axL.text(
        0,
        60.5,
        "s(θ): 상하 슬라이딩, 행정 = 2e = 9 mm",
        color=RED,
        fontsize=11.5,
        ha="center",
        va="bottom",
    )
    axL.annotate(
        "베인 = 별도 부품 (주황)\n1 DOF: y 슬라이딩\n(θ에 종속)",
        xy=(-4.0, 32.0),
        xytext=(-86, 34.0),
        color=COMPRESSION_COLOR,
        fontsize=11.5,
        ha="left",
        va="center",
        linespacing=1.35,
        arrowprops={"arrowstyle": "-", "color": COMPRESSION_COLOR, "lw": 1.2},
    )
    axL.annotate(
        "실린더 벽 관통 슬롯이 베인 안내\n(x이동·회전 구속)",
        xy=(5.0, 40.0),
        xytext=(26, 50.0),
        color=DARK,
        fontsize=11.5,
        ha="left",
        va="center",
        linespacing=1.35,
        arrowprops={"arrowstyle": "-", "color": DARK, "lw": 1.0},
    )
    axL.annotate(
        "팁은 로터 기준선을 따라\n오르내림 (θ=0에서\nR로부터 a_v = 9 mm)",
        xy=(-30.0, 18.0),
        xytext=(-86, 16.0),
        color=COMPRESSION_COLOR,
        fontsize=11.5,
        ha="left",
        va="center",
        linespacing=1.35,
        arrowprops={"arrowstyle": "-", "color": COMPRESSION_COLOR, "lw": 1.2},
    )
    axL.annotate(
        "실린더(빗금 벽): 고정",
        xy=(-34.0, -18.0),
        xytext=(-86, -36.0),
        color=DARK,
        fontsize=11.5,
        ha="left",
        va="center",
        arrowprops={"arrowstyle": "-", "color": DARK, "lw": 1.0},
    )
    axL.annotate(
        "팁이 컷아웃 원 아래에서 로터와 겹침\n(현재 GUI의 투명 겹침 그대로, CAD 확인 대기)",
        xy=(11.0, 6.0),
        xytext=(6, -44.0),
        color=RED,
        fontsize=11.5,
        ha="center",
        va="center",
        linespacing=1.35,
        arrowprops={"arrowstyle": "-", "color": RED, "lw": 1.0},
    )

    # ---------------- Panel B: cylinder-integral fixed vane ----------------
    frame(
        axR,
        "B. 말씀하신 구조로 이해한 것 — 실린더 일체 고정 베인\n"
        "(로터의 구멍·개구부가 고정 베인 둘레를 미끄러짐)",
    )
    draw_body(axR)
    # hatched vane, integral with the wall, hanging down to the drawn tip (y = 25)
    axR.add_patch(
        Rectangle(
            (-half_w, 25.0),
            2 * half_w,
            42.0 - 25.0,
            facecolor=HATCH_BG,
            edgecolor="black",
            linewidth=2.0,
            hatch="///",
            zorder=8,
        )
    )
    # red dashed line from the drawn tip down to y = 13.5 (the a_v result) + X
    axR.plot([0, 0], [25.0, 15.0], color=RED, linestyle=(0, (4, 3)), linewidth=2.2, zorder=9)
    axR.plot(0, 13.5, marker="x", color=RED, markersize=13, markeredgewidth=3.0, zorder=10)
    # rotor-hole sliding stroke arrow (9 mm) to the right of the vane
    axR.annotate(
        "",
        xy=(14.0, 35.0),
        xytext=(14.0, 26.0),
        arrowprops={"arrowstyle": "<->", "color": RED, "lw": 2.2},
        zorder=10,
    )
    axR.annotate(
        "베인 = 실린더와 일체\n(같은 빗금)\n0 DOF, 길이 고정",
        xy=(2.0, 40.0),
        xytext=(-86, 52.0),
        color=DARK,
        fontsize=11.5,
        ha="left",
        va="center",
        linespacing=1.35,
        arrowprops={"arrowstyle": "-", "color": DARK, "lw": 1.0},
    )
    axR.annotate(
        "로터 구멍(Ø16)이 고정 베인 팁 둘레를\n상하 9 mm 슬라이딩 + 최대 ±10.4° 기울어짐",
        xy=(8.0, 31.0),
        xytext=(20, 55.0),
        color=GUIDE_COLOR,
        fontsize=11.5,
        ha="left",
        va="center",
        linespacing=1.35,
        arrowprops={"arrowstyle": "-", "color": GUIDE_COLOR, "lw": 1.0},
    )
    axR.annotate(
        "팁이 항상 구멍 안에\n있으려면 팁 y는 약\n21.5~28.5 mm 필요\n(그림은 y = 25 mm)",
        xy=(-8.0, 27.0),
        xytext=(-86, 12.0),
        color=DARK,
        fontsize=11.5,
        ha="left",
        va="center",
        linespacing=1.35,
        arrowprops={"arrowstyle": "-", "color": DARK, "lw": 1.0},
    )
    axR.annotate(
        "공급 치수 a_v = 9 mm를 그대로 쓰면 팁(θ=0)은\n"
        "y = 13.5 → 로터 몸체 내부(X)\n→ 일체형 해석과 모순, CAD 확인 필요",
        xy=(2.0, 12.0),
        xytext=(4, -44.0),
        color=RED,
        fontsize=11.5,
        ha="center",
        va="center",
        linespacing=1.35,
        arrowprops={"arrowstyle": "-", "color": RED, "lw": 1.0},
    )

    # ---------------- Figure-level title and conclusion ----------------
    fig.suptitle(
        "베인 구조 두 해석 비교 (θ = 0, 상사점) — 고정: 빗금 / 이동: 민무늬", fontsize=16, y=0.965
    )
    fig.text(
        0.5,
        0.035,
        "자유도 결론은 두 모델 모두 동일: 순 자유도 = 크랭크각 θ 1개.  "
        "차이점: A는 베인이 별도 물체(베인-슬롯 마찰·반력, 팁 상하왕복), "
        "B는 베인이 지반의 일부(미끄럼·마찰이 로터 구멍-베인 접촉면에서 발생).",
        fontsize=12.5,
        ha="center",
        va="bottom",
        color=DARK,
    )

    fig.subplots_adjust(left=0.02, right=0.98, top=0.86, bottom=0.09, wspace=0.06)
    fig.savefig(path)
    plt.close(fig)


def render_geometry_master(geometry: RotaryGeometry, path: Path) -> None:
    """Five-panel master plate of the confirmed rotary-compressor geometry.

    Panels: (1) top view at theta=0, (2) axial section on the vane centre
    plane, (3) inlet-lip G1 chain, (4) ledge/recess axial clearance, and
    (5) one swing-bush piece detail, over two summary lines of key numbers.
    Panels (1) and (3) are drawn from the live ``rotor_contour`` so the mouth
    is whatever the profile currently builds.
    """

    RED = "#c0392b"
    GREEN = "#2e7d32"
    BLUE = "#1f5fa8"
    SEAL = "#a9d18e"
    HATCH_FC = "#f2f2f2"
    THIN = "#c8c8c8"

    bore_mm = geometry.cylinder_radius_m / MM
    rotor_mm = geometry.rotor_radius_m / MM
    state = prescribed_state(geometry, 0.0)
    contour = rotor_contour(geometry, 0.0)
    bush = SwingBush()
    bands = AxialBands()
    H = geometry.cylinder_height_m / MM
    ledge = bands.ledge_thickness_m / MM  # 2.4
    ecc = geometry.eccentricity_m / MM  # 4.5
    thin = 14.7  # rotor channel thickness
    recess = 0.5 * (H - thin)  # 3.15 per side

    fig = plt.figure(figsize=(18.4, 15.4), dpi=120)
    fig.suptitle(
        "mochi 로터리 압축기 — 확정 형상 종합도 (단위 mm · θ = 0 상사점 · PHYSICS.md §3, §3.3)",
        fontsize=16,
    )

    def dim_v(ax, x, y0, y1, text, color, dx=1.0, fs=10, ha="left", va="center", rot=0):
        ax.annotate(
            "", xy=(x, y1), xytext=(x, y0), arrowprops=dict(arrowstyle="<->", color=color, lw=1.1)
        )
        ax.text(x + dx, 0.5 * (y0 + y1), text, color=color, fontsize=fs, ha=ha, va=va, rotation=rot)

    def dim_h(ax, y, x0, x1, text, color, dy=1.0, fs=10, va="bottom"):
        ax.annotate(
            "", xy=(x1, y), xytext=(x0, y), arrowprops=dict(arrowstyle="<->", color=color, lw=1.1)
        )
        ax.text(0.5 * (x0 + x1), y + dy, text, color=color, fontsize=fs, ha="center", va=va)

    # ------------------------------------------------------------------
    # (1) Top view at theta = 0
    # ------------------------------------------------------------------
    ax1 = fig.add_axes([0.015, 0.545, 0.30, 0.375])
    ax1.set_aspect("equal")
    ax1.set_xlim(-bore_mm - 20, bore_mm + 14)
    ax1.set_ylim(-bore_mm - 16, bore_mm + 18)
    ax1.axis("off")
    ax1.set_title("① 상면도 (θ = 0, 상사점)", fontsize=13, loc="center")

    ax1.fill(
        [p[0] / MM for p in contour.material],
        [p[1] / MM for p in contour.material],
        facecolor=ROTOR_COLOR,
        edgecolor="none",
        zorder=2,
    )
    groove_x, groove_y = (c / MM for c in state.cutout_center_m)
    for side in (1.0, -1.0):
        xs, ys = _bush_outline_mm(
            state.cutout_center_m[0] + side * bush.piece_shift_m, state.cutout_center_m[1], side
        )
        ax1.fill(xs, ys, facecolor=BUSH_COLOR, edgecolor="black", linewidth=0.8, zorder=3)
    vane_x, vane_y = _vane_outline_mm(geometry, state.vane_tip_m[1])
    ax1.fill(vane_x, vane_y, facecolor=VANE_COLOR, edgecolor="black", linewidth=1.1, zorder=4)
    ax1.add_patch(
        plt.Circle(
            (0.0, 0.0), bore_mm, facecolor="none", edgecolor="black", linewidth=2.2, zorder=5
        )
    )
    for edge in (contour.od_arc, contour.inlet_flat, contour.mouth_path):
        ax1.plot(
            [p[0] / MM for p in edge],
            [p[1] / MM for p in edge],
            color="black",
            linewidth=2.2,
            zorder=5,
        )

    # Dotted vane channel (8*sqrt(2) wide) reaching down into the rotor.
    ch_half = 0.5 * geometry.vane_width_m * sqrt(2.0) / MM
    ax1.add_patch(
        Rectangle(
            (-ch_half, groove_y - 20.0),
            2 * ch_half,
            20.0,
            facecolor="none",
            edgecolor=GUIDE_COLOR,
            linewidth=1.0,
            linestyle=(0, (1, 2)),
            zorder=3,
        )
    )
    for sx in (-ch_half, ch_half):
        ax1.plot(
            [sx, 0.0],
            [groove_y, groove_y - 20.0],
            color=RED,
            linestyle=(0, (2, 3)),
            linewidth=0.9,
            zorder=3,
        )

    # Centre marks.
    rcx, rcy = state.rotor_center_m[0] / MM, state.rotor_center_m[1] / MM
    ax1.plot(rcx, rcy, marker="+", color=RED, markersize=11, markeredgewidth=2, zorder=7)
    ax1.plot(0.0, 0.0, marker="+", color=GUIDE_COLOR, markersize=10, markeredgewidth=2, zorder=7)
    ax1.text(1.4, -3.0, "C", color=GUIDE_COLOR, fontsize=11, fontweight="bold", zorder=7)

    # Ports at +/-30 deg.
    _draw_port_arc(
        ax1, geometry, -30.0, -8.0, "OUT", color="black", label_color="black", label_offset_mm=9.0
    )
    _draw_port_arc(
        ax1, geometry, 8.0, 30.0, "IN", color="black", label_color="black", label_offset_mm=9.0
    )
    for a, ha in ((-30.0, "right"), (30.0, "left")):
        p = port_position(geometry.cylinder_radius_m + 3.0 * MM, a)
        ax1.text(p[0] / MM, p[1] / MM, "30°", color=GUIDE_COLOR, fontsize=9, ha=ha, va="bottom")

    # Vane-width dimension over the top.
    top_y = bore_mm + 3.0
    dim_h(ax1, top_y, -4.0, 4.0, "8 (베인 폭)", GUIDE_COLOR, dy=0.8, fs=9)
    # e = 4.5 near the centre.
    dim_v(ax1, -13.0, 0.0, ecc, "e = 4.5", GUIDE_COLOR, dx=-1.0, fs=9, ha="right")
    # L = 25 from rotor centre up to the groove centre.
    ax1.annotate(
        "",
        xy=(-bore_mm - 6.0, groove_y),
        xytext=(-bore_mm - 6.0, rcy),
        arrowprops=dict(arrowstyle="<->", color=GUIDE_COLOR, lw=1.1),
    )
    ax1.text(
        -bore_mm - 7.0,
        0.5 * (groove_y + rcy),
        "L = 25",
        color=GUIDE_COLOR,
        fontsize=9,
        ha="right",
        va="center",
    )

    ax1.annotate(
        "Ø16 원형 홈",
        xy=(groove_x + 6.0, groove_y + 2.0),
        xytext=(bore_mm * 0.30, bore_mm + 8.0),
        fontsize=9,
        color=GUIDE_COLOR,
        arrowprops=dict(arrowstyle="->", color=GUIDE_COLOR, lw=1.0),
    )
    ax1.annotate(
        "스윙 부시 2조각 → 상세 ⑤",
        xy=(groove_x + bush.piece_outer_radius_m / MM, groove_y),
        xytext=(bore_mm * 0.34, bore_mm + 3.5),
        fontsize=9,
        color=BLUE,
        arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.0),
    )
    ax1.annotate(
        "박육 14.7 채널",
        xy=(-ch_half, groove_y - 14.0),
        xytext=(-bore_mm - 19.0, -6.0),
        fontsize=9,
        color=GUIDE_COLOR,
        arrowprops=dict(arrowstyle="->", color=GUIDE_COLOR, lw=1.0),
    )
    ax1.annotate(
        "8√2 = 11.31",
        xy=(0.0, groove_y - 18.0),
        xytext=(-2.0, groove_y - 27.0),
        fontsize=9,
        color=GUIDE_COLOR,
    )
    ax1.annotate(
        "베인 전체두께 / 턱 구간\n(깊이 15.4 · 25.0 → 상세 ②)",
        xy=(ch_half, groove_y - 8.0),
        xytext=(ch_half + 6.0, groove_y - 16.0),
        fontsize=9,
        color=GUIDE_COLOR,
        arrowprops=dict(arrowstyle="->", color=GUIDE_COLOR, lw=1.0),
    )
    ax1.annotate(
        "Ø77 실린더 보어",
        xy=(-bore_mm * 0.72, -bore_mm * 0.72),
        xytext=(-bore_mm - 6.0, -bore_mm - 6.0),
        fontsize=9,
        color=BLUE,
        arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.0),
    )
    ax1.annotate(
        "Ø68 로터 OD",
        xy=(-rotor_mm * 0.66, -rotor_mm * 0.66),
        xytext=(-bore_mm - 6.0, -bore_mm - 10.5),
        fontsize=9,
        color=BLUE,
        arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.0),
    )
    ax1.text(
        bore_mm * 0.02,
        bore_mm + 12.5,
        "흡입측만 OD 평면부 (34.8°→13.4°)\n토출측은 OD 접선 13° → 상세 ③",
        color=RED,
        fontsize=8.5,
        ha="left",
        va="bottom",
    )

    # ------------------------------------------------------------------
    # (2) Axial section on the vane centre plane
    # ------------------------------------------------------------------
    ax2 = fig.add_axes([0.35, 0.545, 0.635, 0.375])
    ax2.set_aspect("equal")
    ax2.set_xlim(-6, 74)
    ax2.set_ylim(-16, 40)
    ax2.axis("off")
    ax2.set_title("② 측면 단면 (베인 중심면 x = 0)", fontsize=13)

    plate_t = 9.0
    # End plates (top and bottom), hatched.
    ax2.add_patch(
        Rectangle(
            (0, H), 62, plate_t, facecolor=HATCH_FC, edgecolor="black", hatch="////", linewidth=1.0
        )
    )
    ax2.add_patch(
        Rectangle(
            (0, -plate_t),
            62,
            plate_t,
            facecolor=HATCH_FC,
            edgecolor="black",
            hatch="////",
            linewidth=1.0,
        )
    )
    # Rotor solid block (full height).
    ax2.add_patch(Rectangle((6, 0), 26, H, facecolor=ROTOR_COLOR, edgecolor="black", linewidth=1.2))
    ax2.text(19, H / 2, "로터 (전체 두께 21)", ha="center", va="center", fontsize=10)
    # Thinned rotor channel (박육) toward the groove.
    ax2.add_patch(
        Rectangle((38, recess), 12, thin, facecolor=THIN, edgecolor="black", linewidth=1.0)
    )
    ax2.text(44, H / 2, "박육", ha="center", va="center", fontsize=9)
    # Green sealed bands (4 MPa) top and bottom of the ledge region.
    for yb in (H - recess, 0.0):
        ax2.add_patch(Rectangle((38, yb), 20, recess, facecolor=SEAL, edgecolor="none", zorder=1))
    # Vane full block down into the groove (right of channel), hatched wall.
    ax2.add_patch(
        Rectangle(
            (50, recess),
            12,
            thin,
            facecolor=HATCH_FC,
            edgecolor="black",
            hatch="////",
            linewidth=1.0,
        )
    )

    # Dimensions.
    dim_v(ax2, 2.0, 0.0, H, "21", GUIDE_COLOR, dx=-1.0, ha="right")
    dim_v(
        ax2, 35.5, recess, H - recess, "14.7", BLUE, dx=-0.8, fs=9, rot=90, va="center", ha="right"
    )
    dim_v(ax2, 66.0, ledge, H - ledge, "16.2", BLUE, dx=0.8, fs=9)
    dim_h(ax2, 36.0, 20.0, 62.0, "25.0", GUIDE_COLOR, dy=0.8, fs=10)
    dim_h(ax2, 32.5, 38.0, 62.0, "15.4", GUIDE_COLOR, dy=0.8, fs=10)

    # detail (4) callout.
    ax2.add_patch(
        plt.Circle(
            (49, H - 1.0),
            5.0,
            facecolor="none",
            edgecolor=RED,
            linewidth=1.2,
            linestyle=(0, (4, 3)),
            zorder=6,
        )
    )
    ax2.annotate(
        "상세 ④",
        xy=(52, H + 1.5),
        xytext=(24, 34.0),
        fontsize=10,
        color=RED,
        arrowprops=dict(arrowstyle="->", color=RED, lw=1.1),
    )
    # a_v = 9 with R -> tip along the bottom.
    ax2.annotate(
        "R",
        xy=(40, -plate_t - 1.0),
        xytext=(40, -plate_t - 3.0),
        color=RED,
        fontsize=10,
        ha="center",
        va="top",
    )
    dim_h(ax2, -plate_t - 5.0, 40.0, 58.0, "a_v = 9", RED, dy=-3.0, fs=11, va="top")

    # Green legend under panel (2).
    fig.text(0.395, 0.520, "■", color=SEAL, fontsize=15, ha="left", va="center")
    fig.text(
        0.415,
        0.520,
        "4 MPa(절대) 밀폐 공간 — 로터 21부·부시·베인이 폐쇄",
        fontsize=10,
        ha="left",
        va="center",
    )

    # ------------------------------------------------------------------
    # (3) Inlet-lip G1 chain
    # ------------------------------------------------------------------
    ax3 = fig.add_axes([0.015, 0.10, 0.31, 0.35])
    ax3.set_aspect("equal")
    ax3.axis("off")
    ax3.set_title("③ 흡입측 립 (OD 평면부 + G1 체인)", fontsize=13)

    mouth = [(x / MM, y / MM) for x, y in contour.mouth_path]
    inlet_lip_len = 1 + 41 + 1 + 31
    lip = mouth[: inlet_lip_len + 2]
    gcx, gcy = groove_x, groove_y
    arc = mouth[1:42]
    straight = mouth[41:43]
    blend = mouth[42:inlet_lip_len]
    inlet_touch = mouth[inlet_lip_len]
    inlet_angle = degrees(atan2(inlet_touch[0] - gcx, inlet_touch[1] - gcy))

    for edge in (contour.od_arc, contour.inlet_flat):
        ax3.plot([p[0] / MM for p in edge], [p[1] / MM for p in edge], color="black", linewidth=2.2)
    ax3.plot([p[0] for p in arc], [p[1] for p in arc], color=RED, linewidth=2.6)
    ax3.plot([p[0] for p in straight], [p[1] for p in straight], color="black", linewidth=2.6)
    ax3.plot([p[0] for p in blend], [p[1] for p in blend], color=GREEN, linewidth=2.6)
    # Rest of the mouth (down the groove and back) faint for context.
    ax3.plot(
        [p[0] for p in mouth[inlet_lip_len - 1 :]],
        [p[1] for p in mouth[inlet_lip_len - 1 :]],
        color="black",
        linewidth=1.4,
        alpha=0.55,
    )

    # Dotted groove-centre ray through the inlet touch.
    ray_dx, ray_dy = inlet_touch[0] - gcx, inlet_touch[1] - gcy
    ax3.plot(
        [gcx - 0.4 * ray_dx, inlet_touch[0] + 1.6 * ray_dx],
        [gcy - 0.4 * ray_dy, inlet_touch[1] + 1.6 * ray_dy],
        color=RED,
        linestyle=(0, (2, 3)),
        linewidth=1.0,
    )

    ax3.plot(*mouth[0], "o", color=BLUE, markersize=7)
    ax3.plot(*inlet_touch, "o", color=BLUE, markersize=7)
    ax3.plot(*mouth[41], "o", color=RED, markersize=6)
    ax3.plot(*mouth[42], "o", color=GREEN, markersize=6)
    ax3.plot(*inlet_touch, "o", color=GREEN, markersize=6)

    x0 = min(p[0] for p in lip) - 8.0
    x1 = max(p[0] for p in lip) + 6.0
    y0 = min(p[1] for p in lip) - 6.0
    y1 = max(p[1] for p in lip) + 5.5
    ax3.set_xlim(x0, x1)
    ax3.set_ylim(y0, y1)

    ax3.annotate(
        "B: 평면 끝 13.4° (R33.657)",
        xy=mouth[0],
        xytext=(mouth[0][0] + 1.4, mouth[0][1] + 3.2),
        fontsize=9,
        color=BLUE,
        ha="left",
        arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.0),
    )
    ax3.annotate(
        "R1.5 × 102.2°",
        xy=arc[18],
        xytext=(arc[18][0] - 4.0, arc[18][1] + 4.2),
        fontsize=9,
        color=RED,
        ha="right",
        arrowprops=dict(arrowstyle="->", color=RED, lw=1.0),
    )
    ax3.annotate(
        "직선 0.399",
        xy=straight[0],
        xytext=(straight[0][0] - 6.5, straight[0][1] - 1.0),
        fontsize=9,
        color="black",
        ha="right",
        arrowprops=dict(arrowstyle="->", color="black", lw=1.0),
    )
    ax3.annotate(
        "R1.0 × 52.62°",
        xy=blend[14],
        xytext=(blend[14][0] + 3.0, blend[14][1] - 1.2),
        fontsize=9,
        color=GREEN,
        arrowprops=dict(arrowstyle="->", color=GREEN, lw=1.0),
    )
    ax3.annotate(
        f"D: 홈 블렌드 {abs(inlet_angle):.2f}°",
        xy=inlet_touch,
        xytext=(inlet_touch[0] - 5.5, inlet_touch[1] - 4.0),
        fontsize=9,
        color=GREEN,
        ha="right",
        arrowprops=dict(arrowstyle="->", color=GREEN, lw=1.0),
    )

    # ------------------------------------------------------------------
    # (4) Ledge / recess axial clearance
    # ------------------------------------------------------------------
    ax4 = fig.add_axes([0.355, 0.10, 0.29, 0.35])
    ax4.set_aspect("equal")
    ax4.set_xlim(-1, 25)
    ax4.set_ylim(-1, 13)
    ax4.axis("off")
    ax4.set_title("④ 턱-리세스 유격 (축방향 상세)", fontsize=13)

    # Top end plate with a 2.4 ledge (턱) stepping down over the vane.
    ax4.add_patch(
        Rectangle(
            (0, 9.0), 20, 3.5, facecolor=HATCH_FC, edgecolor="black", hatch="////", linewidth=1.0
        )
    )
    ax4.add_patch(
        Rectangle(
            (8, 6.6), 12, 2.4, facecolor=HATCH_FC, edgecolor="black", hatch="////", linewidth=1.0
        )
    )  # 2.4 ledge
    # Green sealed L-band over the recess (3.15 deep, above the rotor thin part).
    ax4.add_patch(Rectangle((0, 3.45), 20, 3.15, facecolor=SEAL, edgecolor="none"))
    ax4.add_patch(Rectangle((0, 3.45), 8, 5.55, facecolor=SEAL, edgecolor="none"))
    # Rotor thin part (박육부) below.
    ax4.add_patch(Rectangle((0, 0.0), 12, 3.45, facecolor=THIN, edgecolor="black", linewidth=1.0))
    ax4.text(6, 1.6, "로터 박육부", ha="center", va="center", fontsize=9)
    # Vane block (hatched) at right.
    ax4.add_patch(
        Rectangle(
            (16, 0.0), 4, 6.6, facecolor=HATCH_FC, edgecolor="black", hatch="////", linewidth=1.0
        )
    )
    ax4.text(18, 3.0, "베인", ha="center", va="center", fontsize=8, rotation=90)

    dim_v(ax4, 1.4, 3.45, 6.6, "3.15", BLUE, dx=0.4, fs=9)
    ax4.annotate(
        "",
        xy=(21.0, 9.0),
        xytext=(21.0, 6.6),
        arrowprops=dict(arrowstyle="<->", color=BLUE, lw=1.1),
    )
    ax4.text(21.4, 7.8, "2.4 턱", color=BLUE, fontsize=9, ha="left", va="center")
    dim_v(ax4, 10.5, 3.45, 4.2, "0.75", BLUE, dx=0.4, fs=9)

    # ------------------------------------------------------------------
    # (5) Swing-bush piece detail
    # ------------------------------------------------------------------
    ax5 = fig.add_axes([0.665, 0.10, 0.31, 0.35])
    ax5.set_aspect("equal")
    ax5.axis("off")
    ax5.set_title("⑤ 스윙 부시 (곡면 106.1° + R0.5 블렌드)", fontsize=13)

    r_out = bush.piece_outer_radius_m / MM
    # Two bush pieces hugging a central vane block, drawn large.
    vane_hw = 0.5 * geometry.vane_width_m / MM
    for side in (1.0, -1.0):
        xs, ys = _bush_outline_mm(side * bush.piece_shift_m, 0.0, side)
        ax5.fill(xs, ys, facecolor="#a9b8d0", edgecolor="black", linewidth=1.3, zorder=2)
    ax5.add_patch(
        Rectangle(
            (-vane_hw, -r_out - 3.5),
            2 * vane_hw,
            2 * r_out + 7.0,
            facecolor=VANE_COLOR,
            edgecolor="black",
            linewidth=1.2,
            zorder=3,
        )
    )
    ax5.plot(0.0, 0.0, marker="+", color=RED, markersize=13, markeredgewidth=2.4, zorder=4)
    ax5.text(0.9, 0.4, "H", color=RED, fontsize=11, fontweight="bold", zorder=4)

    ax5.set_xlim(-r_out - 5, r_out + 5)
    ax5.set_ylim(-r_out - 5, r_out + 5)

    # Curved-face arc angle marker (106.1 deg) on the right piece.
    half_arc = degrees(bush.half_arc_rad())
    aa = [radians(90 - a) for a in [-(half_arc) + 2 * half_arc * i / 40 for i in range(41)]]
    ax5.plot([2.4 * cos(t) for t in aa], [2.4 * sin(t) for t in aa], color=RED, linewidth=1.0)

    ax5.annotate(
        "곡면 106.1°",
        xy=(r_out * 0.7, r_out * 0.55),
        xytext=(r_out * 0.55, r_out + 3.0),
        fontsize=9,
        color=RED,
        ha="left",
        arrowprops=dict(arrowstyle="->", color=RED, lw=1.0),
    )
    ax5.annotate(
        "R0.5 블렌드",
        xy=(-vane_hw - 0.4, r_out - 1.6),
        xytext=(-r_out - 3.5, r_out + 2.0),
        fontsize=9,
        color=GREEN,
        arrowprops=dict(arrowstyle="->", color=GREEN, lw=1.0),
    )
    ax5.annotate(
        "평면 3.990\n(접촉 길이 11.94)",
        xy=(-vane_hw, -r_out * 0.55),
        xytext=(-r_out - 4.5, -r_out - 2.0),
        fontsize=9,
        color=BLUE,
        arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.0),
    )
    ax5.annotate(
        "곡면 R7.970\n(유막 0.010)",
        xy=(r_out * 0.7, -r_out * 0.6),
        xytext=(r_out - 1.0, -r_out - 2.0),
        fontsize=9,
        color=BLUE,
        arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.0),
    )

    # ------------------------------------------------------------------
    # Bottom summary lines
    # ------------------------------------------------------------------
    fig.text(
        0.5,
        0.055,
        "실린더 Ø77 · 로터 Ø68 · 편심 4.5 · 원형 홈 Ø16 (L = 25) · "
        "베인 폭 8 (전체두께 보어~15.4, 팁 25.0, 턱 2.4×2 → 빈틈 16.2) · "
        "로터 박육 14.7 (리세스 3.15/측, 채널 폭 8√2 = 11.31) · 축 높이 21 · 포트 ±30°",
        ha="center",
        va="center",
        fontsize=11,
    )
    fig.text(
        0.5,
        0.022,
        "스윙 부시: 조각 OD R7.970 (Ø15.94, 곡면 106.1°) · 평면부 3.990 (접촉 11.94) · "
        "사이는 R0.5 블렌드 → 조각 중심을 ±0.020 이동 → 평면-베인 0.010 / 곡면-홈 0.010 "
        "(면내 유격 0.020의 중앙) · 축방향 0.0085/편측 → 높이 20.983",
        ha="center",
        va="center",
        fontsize=11,
        color=BLUE,
    )

    fig.savefig(path)
    plt.close(fig)


def render_chamber_case1_seal_over(geometry: RotaryGeometry, path: Path) -> None:
    """Case 1 seal-over near TDC: the rotor-cylinder contact point hides under
    the vane and the suction/discharge chambers merge into one connected
    (purple) region (circular-rotor approximation, PHYSICS section 3.1/3.2)."""

    bore = geometry.cylinder_radius_m / MM
    vane_w = geometry.vane_width_m / MM
    half_w = 0.5 * vane_w
    half_deg = degrees(seal_over_half_angle_rad(geometry))
    theta_r = radians(3.0)

    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(14.0, 7.7), dpi=150)
    fig.suptitle(
        "케이스 1 — 밀봉 전환(seal-over) 구간: "
        "접촉점이 베인 아래로 들어가면 두 영역이 하나로 합쳐진다",
        fontsize=14,
    )

    # -- left: the normal (separated) case, two chambers split at contact T ---
    _draw_dead_centre_panel(ax_l, geometry, 180.0, "정상 분리 — θ = 180°", merged=False)
    ax_l.set_xlim(-bore - 13, bore + 13)
    ax_l.set_ylim(-bore - 22, bore + 12)
    ax_l.annotate(
        "벽 1: 베인",
        xy=(1.0, 29.0),
        xytext=(-30.0, 47.0),
        fontsize=11,
        color="black",
        ha="left",
        arrowprops={"arrowstyle": "->", "color": "black", "lw": 1.2},
    )
    ax_l.text(
        0.0,
        -bore - 6.0,
        "벽 2: 로터-실린더 접촉점 T\n(밀봉점, 여기서 두 영역이 나뉨)",
        fontsize=11,
        color="#c0392b",
        ha="center",
        va="top",
    )
    ax_l.annotate(
        "흡입부 (IN 쪽)",
        xy=(34.0, -16.0),
        xytext=(24.0, -bore - 8.0),
        fontsize=12,
        color=SUCTION_COLOR,
        ha="left",
        arrowprops={"arrowstyle": "->", "color": SUCTION_COLOR, "lw": 1.4},
    )
    ax_l.annotate(
        "토출부 (OUT 쪽)",
        xy=(-33.0, -18.0),
        xytext=(-bore - 12.0, -bore + 2.0),
        fontsize=12,
        color=COMPRESSION_COLOR,
        ha="left",
        arrowprops={"arrowstyle": "->", "color": COMPRESSION_COLOR, "lw": 1.4},
    )

    # -- right: the merged (seal-over) case, one connected purple region -------
    _draw_dead_centre_panel(
        ax_r, geometry, 3.0, "케이스 1: 병합 구간 — θ = 3° (< 5.96°)", merged=True
    )
    ax_r.set_xlim(-bore - 13, bore + 13)
    ax_r.set_ylim(-bore - 22, bore + 12)

    # red dashed +-5.96 deg merge window, from centre C out to the bore
    for sgn in (-1.0, 1.0):
        ax_r.plot(
            [0.0, bore * sin(sgn * radians(half_deg))],
            [0.0, bore * cos(sgn * radians(half_deg))],
            color="#c0392b",
            linewidth=1.1,
            linestyle=(0, (5, 4)),
            zorder=6,
        )

    tx, ty = bore * sin(theta_r), bore * cos(theta_r)
    ax_r.annotate(
        "접촉점 T가 베인 폭(8 mm) 안에 숨음\n→ 벽 2가 사라짐",
        xy=(tx, ty),
        xytext=(-bore - 12.0, bore + 4.0),
        fontsize=11,
        color="#c0392b",
        ha="left",
        arrowprops={"arrowstyle": "->", "color": "#c0392b", "lw": 1.2},
    )
    ax_r.annotate(
        "빨간 점선: ±5.96° 병합 구간",
        xy=(-1.6, 24.0),
        xytext=(-bore - 12.0, 20.0),
        fontsize=11,
        color="#c0392b",
        ha="left",
        arrowprops={"arrowstyle": "->", "color": "#c0392b", "lw": 1.2},
    )
    ax_r.annotate(
        "전체가 하나의 보라색 영역:\n오른쪽(IN)과 왼쪽(OUT)이\n아래로 돌아 그대로 연결됨",
        xy=(26.0, -27.0),
        xytext=(4.0, -bore - 8.0),
        fontsize=11,
        color=MERGED_COLOR,
        ha="left",
        arrowprops={"arrowstyle": "->", "color": MERGED_COLOR, "lw": 1.4},
    )

    # zoom inset: contact point T tucked under the vane width
    axz = ax_r.inset_axes([0.63, 0.70, 0.40, 0.26])
    axz.set_title("확대: T가 베인 아래", fontsize=10)
    zx = [bore * sin(radians(a)) for a in range(-32, 33)]
    zy = [bore * cos(radians(a)) for a in range(-32, 33)]
    vane_top = 37.0
    axz.add_patch(
        Rectangle(
            (-half_w, 30.0),
            vane_w,
            vane_top - 30.0,
            facecolor=VANE_COLOR,
            edgecolor="black",
            linewidth=1.0,
            zorder=3,
        )
    )
    for sgn in (-1.0, 1.0):
        axz.plot(
            [sgn * half_w, sgn * half_w],
            [33.5, vane_top],
            color="white",
            linewidth=1.1,
            linestyle=(0, (2, 2)),
            zorder=4,
        )
    axz.plot(zx, zy, color="black", linewidth=2.0, zorder=2)
    axz.plot(tx, ty, "o", color="#c0392b", markersize=9, zorder=5)
    axz.set_xlim(-14.0, 14.0)
    axz.set_ylim(33.5, 40.2)
    axz.set_xticks([])
    axz.set_yticks([])

    fig.tight_layout(rect=(0.0, 0.03, 1.0, 0.95))
    fig.savefig(path)
    plt.close(fig)


def render_chamber_case2_gap(geometry: RotaryGeometry, path: Path) -> None:
    """Case 2 -- eccentricity below the tangency condition opens a rotor/
    cylinder gap that connects the two chambers (no valid split)."""

    bore = geometry.cylinder_radius_m / MM
    rotor = geometry.rotor_radius_m / MM
    half_deg = degrees(seal_over_half_angle_rad(geometry))
    theta_deg = 90.0
    theta = radians(theta_deg)

    # Reduced eccentricity: a 0.5 mm radial gap remains at the intended contact.
    gap_geom = replace(geometry, eccentricity_m=4.0 * MM)

    red = "#c0392b"
    out_text = "#b0642c"

    def draw_panel(ax, geom, fills, title):
        ecc = geom.eccentricity_m / MM
        ax.set_aspect("equal")
        ax.set_xlim(-bore - 6, bore + 6)
        ax.set_ylim(-bore - 12, bore + 8)
        ax.axis("off")
        ax.set_title(title, fontsize=13, pad=14)

        for a0, a1, color in fills:
            _fill_crescent_sector(ax, geom, theta, a0, a1, color)

        orbit_x, orbit_y = ecc * sin(theta), ecc * cos(theta)
        ax.add_patch(
            plt.Circle(
                (orbit_x, orbit_y),
                rotor,
                facecolor=ROTOR_COLOR,
                edgecolor="black",
                linewidth=1.5,
                zorder=2,
            )
        )
        ax.add_patch(
            plt.Circle(
                (0.0, 0.0), bore, facecolor="none", edgecolor="black", linewidth=2.2, zorder=3
            )
        )

        state = prescribed_state(geom, theta)
        half_width = 0.5 * geom.vane_width_m / MM
        vane_tip_y = state.vane_tip_m[1] / MM
        ax.add_patch(
            Rectangle(
                (-half_width, vane_tip_y),
                2.0 * half_width,
                bore - vane_tip_y,
                facecolor=VANE_COLOR,
                edgecolor="black",
                linewidth=1.0,
                zorder=4,
            )
        )

        # Bore centre C (black) and rotor centre R (guide grey), no text labels.
        ax.plot(0.0, 0.0, marker="+", color="black", markersize=11, markeredgewidth=2.2, zorder=6)
        ax.plot(
            orbit_x,
            orbit_y,
            marker="+",
            color=GUIDE_COLOR,
            markersize=11,
            markeredgewidth=2.2,
            zorder=6,
        )

        for angle_deg, text, ha in ((-40.0, "OUT", "right"), (40.0, "IN", "left")):
            tick_in = port_position(geom.cylinder_radius_m - 1.0 * MM, angle_deg)
            tick_out = port_position(geom.cylinder_radius_m + 2.5 * MM, angle_deg)
            ax.plot(
                [tick_in[0] / MM, tick_out[0] / MM],
                [tick_in[1] / MM, tick_out[1] / MM],
                color="black",
                linewidth=2,
            )
            label = port_position(geom.cylinder_radius_m + 5.0 * MM, angle_deg)
            ax.text(
                label[0] / MM,
                label[1] / MM,
                text,
                fontweight="bold",
                fontsize=12,
                ha=ha,
                va="center",
            )

    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(15.0, 8.1), dpi=150)
    fig.suptitle(
        "케이스 2 — 편심이 접선 조건보다 작으면 로터·실린더 사이 틈으로 두 영역이 연결된다",
        fontsize=15,
    )

    # ---- LEFT panel: sealed / tangent contact ----------------------------
    draw_panel(
        ax_l,
        geometry,
        [(half_deg, theta_deg, SUCTION_COLOR), (theta_deg, 360.0 - half_deg, COMPRESSION_COLOR)],
        "정상 — e = 4.5 mm = (77 - 68)/2, θ = 90°",
    )
    # Contact point T rides the bore at the crank angle (3 o'clock).
    ax_l.plot(bore * sin(theta), bore * cos(theta), "o", color=red, markersize=10, zorder=7)
    ax_l.annotate(
        "흡입부 (IN 쪽)",
        xy=(bore * 0.72, bore * 0.62),
        xytext=(bore + 2, bore * 0.72),
        fontsize=12,
        color=SUCTION_COLOR,
        ha="left",
        va="center",
        arrowprops={"arrowstyle": "->", "color": SUCTION_COLOR, "lw": 1.4},
    )
    ax_l.annotate(
        "토출부 (OUT 쪽)",
        xy=(-bore * 0.80, -bore * 0.45),
        xytext=(-bore - 4, -bore - 6),
        fontsize=12,
        color=out_text,
        ha="left",
        va="center",
        arrowprops={"arrowstyle": "->", "color": out_text, "lw": 1.4},
    )
    ax_l.annotate(
        "접선(밀봉) 접촉점 T",
        xy=(bore * sin(theta), bore * cos(theta)),
        xytext=(bore * 0.10, -bore - 10),
        fontsize=12,
        color=red,
        ha="center",
        va="center",
        arrowprops={"arrowstyle": "->", "color": red, "lw": 1.4},
    )

    # ---- RIGHT panel: unsealed / connected -------------------------------
    draw_panel(
        ax_r,
        gap_geom,
        [(half_deg, 360.0 - half_deg, MERGED_COLOR)],
        "케이스 2: e = 4.0 mm → 0.5 mm 틈, θ = 90°",
    )
    ax_r.annotate(
        "틈으로 위·아래가 연결되어\n전체가 하나의 보라색 영역",
        xy=(-bore * 0.86, -bore * 0.32),
        xytext=(-bore - 4, -bore - 8),
        fontsize=12,
        color=MERGED_COLOR,
        ha="left",
        va="center",
        arrowprops={"arrowstyle": "->", "color": MERGED_COLOR, "lw": 1.4},
    )
    ax_r.annotate(
        "로터가 실린더에 닿지 않음:\n반경 방향 0.5 mm 틈",
        xy=(bore * 0.985, -bore * 0.20),
        xytext=(bore * 0.30, -bore - 8),
        fontsize=12,
        color=red,
        ha="left",
        va="center",
        arrowprops={"arrowstyle": "->", "color": red, "lw": 1.4},
    )

    # ---- Inset zoom of the 0.5 mm gap at the intended contact -------------
    axins = ax_r.inset_axes([0.86, 0.51, 0.19, 0.49])
    a_lo, a_hi = 82.0, 98.0
    steps = 60
    outer_x, outer_y, inner_x, inner_y = [], [], [], []
    ecc_g = gap_geom.eccentricity_m / MM
    for i in range(steps + 1):
        a = radians(a_lo + (a_hi - a_lo) * i / steps)
        outer_x.append(bore * sin(a))
        outer_y.append(bore * cos(a))
        off = a - theta
        rho = ecc_g * cos(off) + sqrt(max(rotor * rotor - ecc_g * ecc_g * sin(off) ** 2, 0.0))
        inner_x.append(rho * sin(a))
        inner_y.append(rho * cos(a))
    x_min = 37.35
    # Rotor interior (grey) to the left of the inner arc.
    axins.fill(
        inner_x + [x_min, x_min],
        inner_y + [inner_y[-1], inner_y[0]],
        facecolor=ROTOR_COLOR,
        edgecolor="none",
        zorder=0,
    )
    # The connecting gap (purple) between rotor and bore.
    axins.fill(
        outer_x + inner_x[::-1],
        outer_y + inner_y[::-1],
        facecolor=MERGED_COLOR,
        edgecolor="none",
        zorder=1,
    )
    axins.plot(inner_x, inner_y, color="black", linewidth=1.6, zorder=2)
    axins.plot(outer_x, outer_y, color="black", linewidth=2.2, zorder=2)
    axins.set_xlim(x_min, 38.9)
    axins.set_ylim(-5.6, 5.6)
    axins.set_xticks([])
    axins.set_yticks([])
    for sp in axins.spines.values():
        sp.set_edgecolor("black")
        sp.set_linewidth(1.2)
    axins.set_title("확대: 틈새", fontsize=11)
    axins.annotate(
        "",
        xy=(38.12, 4.2),
        xytext=(38.12, -4.2),
        arrowprops={"arrowstyle": "<->", "color": red, "lw": 2.2},
        zorder=3,
    )
    axins.text(37.42, 0.0, "0.5 mm", color=red, fontsize=11, ha="left", va="center", zorder=3)

    ax_r.indicate_inset_zoom(axins, edgecolor="#8a8a8a", alpha=0.85, linewidth=1.0)

    fig.tight_layout(rect=(0.0, 0.02, 1.0, 0.95))
    fig.savefig(path)
    plt.close(fig)


# --------------------------------------------------------------------------
# Figure: indicated P-V diagram and indicated power
# --------------------------------------------------------------------------


def render_indicated_pv_diagram(curve: PressureCurve, path: Path) -> None:
    """Pressure-volume indicator loop and indicated power (PHYSICS.md 3.5)."""

    geometry = curve.geometry
    trace = curve.trace
    work = indicated_work_j(geometry, trace=trace)

    seal_entry = radians(curve.seal_over_entry_deg)
    steps = 400
    angles = [seal_entry * i / steps for i in range(steps + 1)]
    compression_v = [trace.compression_at(a) * 1.0e6 for a in angles]
    suction_v = [trace.suction_at(a) * 1.0e6 for a in angles]
    compression_p = [
        port_timed_pressures(geometry, a, trace=trace).compression_pressure_pa * PA_TO_MPA
        for a in angles
    ]
    suction_p = [
        port_timed_pressures(geometry, a, trace=trace).suction_pressure_pa * PA_TO_MPA
        for a in angles
    ]
    valved = valved_cycle(geometry, trace=trace, samples=2880)
    valved_p = [valved.compression_at(a) * PA_TO_MPA for a in angles]

    # Closed loop: the suction branch (clearance -> max volume) followed by the
    # compression branch (max -> clearance volume) forms the indicator polygon.
    fig, ax = plt.subplots(figsize=(9.6, 6.6), dpi=140)
    ax.fill(
        suction_v + compression_v,
        suction_p + compression_p,
        facecolor=COMPRESSION_COLOR,
        alpha=0.15,
        edgecolor="none",
        zorder=1,
    )
    ax.fill_between(
        compression_v,
        compression_p,
        valved_p,
        facecolor="#c0392b",
        alpha=0.18,
        edgecolor="none",
        zorder=2,
    )
    ax.plot(
        compression_v,
        compression_p,
        color=COMPRESSION_COLOR,
        linewidth=2.2,
        label="압축/토출 (이상)",
    )
    ax.plot(
        compression_v,
        valved_p,
        color="#c0392b",
        linewidth=1.6,
        linestyle="--",
        label="리드밸브 과압 (§3.8, 별도 성능 항)",
        zorder=3,
    )
    ax.plot(suction_v, suction_p, color=SUCTION_COLOR, linewidth=2.2, label="흡입 챔버")
    for pressure in (curve.suction_port_mpa, curve.discharge_port_mpa):
        ax.axhline(pressure, color="#9aa0a8", linestyle=":", linewidth=1.0)

    ax.set_xlabel("챔버 체적 V (cm³)")
    ax.set_ylabel("챔버 압력 p (MPa, 절대압)")
    ax.set_title("지시 선도 (P-V) — 포트-타이밍 실형상 체적 규칙 (PHYSICS.md §3.5)", color="none")
    ax.grid(True, color="#e2e5ea", linewidth=0.8)
    _legend_above(ax, ncol=3, fontsize=8.5)
    ax.set_xlim(-1.0, max(compression_v) * 1.05)
    ax.set_ylim(0.0, curve.residual_mpa * 1.05)
    ax.text(
        0.03,
        0.97,
        f"지시 일 W = {work.net_work_j:.1f} J → 지시동력 {work.power_w:.0f} W (등엔트로피-이상)\n"
        f"리드밸브 과압 (§3.8, 별도·loads 미반영): "
        f"+{curve.overpressure_power_w:.0f} → {curve.valve_indicated_power_w:.0f} W (붉은 음영)\n"
        f"= 압축 {work.compression_work_j:.1f} + 흡입 {work.suction_work_j:.1f} J\n"
        "실측은 여기에 흡입 스로틀·열전달 추가로 +15~40%",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9.5,
        bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "#c7ccd4"},
    )
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# --------------------------------------------------------------------------
# Figure: swing-bush lubrication film
# --------------------------------------------------------------------------


def render_bush_film_pressure(curve: PressureCurve, path: Path) -> None:
    """Swing-bush film pressure distribution (PHYSICS.md 3.6)."""

    geometry = curve.geometry
    theta_rep_deg = 200.0
    state = film_state(geometry, radians(theta_rep_deg), trace=curve.trace)

    fig, ax = plt.subplots(figsize=(7.8, 5.2), dpi=140)
    faces = (
        ("IN 평면", state.in_piece.flat, SUCTION_COLOR, "-"),
        ("IN 곡면", state.in_piece.curved, SUCTION_COLOR, "--"),
        ("OUT 평면", state.out_piece.flat, COMPRESSION_COLOR, "-"),
        ("OUT 곡면", state.out_piece.curved, COMPRESSION_COLOR, "--"),
    )
    for label, face, color, style in faces:
        ax.plot(
            [0.0, face.length_m * 1.0e3],
            [face.inlet_pressure_pa * PA_TO_MPA, face.outlet_pressure_pa * PA_TO_MPA],
            color=color,
            linestyle=style,
            linewidth=2.0,
            label=label,
        )
    ax.set_xlabel("막 위치 x (mm, 챔버단 → 리세스단)")
    ax.set_ylabel("막 압력 p (MPa)")
    ax.set_title(
        f"스윙 부시 막 압력 분포 (θ = {theta_rep_deg:.0f}°) — 리세스 4 MPa 선형 (§3.6)",
        fontsize=10,
        color="none",
    )
    ax.grid(True, color="#e2e5ea", linewidth=0.7)
    _legend_above(ax, ncol=4, fontsize=8)
    _save(fig, path)


def render_bush_film_velocity(curve: PressureCurve, path: Path) -> None:
    """Swing-bush face sliding velocities over a revolution (PHYSICS.md 3.6)."""

    geometry = curve.geometry
    degrees_axis = list(range(0, 361, 2))
    flat_u = [flat_slide_velocity(geometry, radians(d)) for d in degrees_axis]
    curved_u = [curved_slide_velocity(geometry, radians(d)) for d in degrees_axis]

    fig, ax = plt.subplots(figsize=(7.8, 5.2), dpi=140)
    ax.plot(degrees_axis, flat_u, color="#2e7d32", linewidth=1.8, label="평면(병진)")
    ax.plot(degrees_axis, curved_u, color="#7a4fb0", linewidth=1.8, label="곡면(진동)")
    ax.axhline(0.0, color="#c7ccd4", linewidth=0.8)
    ax.set_xlim(0, 360)
    ax.set_xticks(range(0, 361, 90))
    ax.set_xlabel("크랭크 각 θ (deg)")
    ax.set_ylabel("활주속도 U (m/s)")
    ax.set_title("면별 활주속도 U(θ) — 부시 병진·로터 진동 (§3.6)", fontsize=10, color="none")
    ax.grid(True, color="#e2e5ea", linewidth=0.7)
    _legend_above(ax, ncol=2)
    _save(fig, path)


def render_bush_film_shear(curve: PressureCurve, path: Path) -> None:
    """Swing-bush face shear traction with the cycle friction summary (PHYSICS.md 3.6)."""

    geometry = curve.geometry
    trace = curve.trace
    mu = LUBRICANT_VISCOSITY_PA_S
    state = film_state(geometry, radians(200.0), trace=trace)
    film_h = state.in_piece.flat.film_thickness_m
    friction_w = friction_power_cycle_w(geometry, samples=180, trace=trace)
    work = indicated_work_j(geometry, trace=trace)

    degrees_axis = list(range(0, 361, 2))
    flat_tau = [mu * flat_slide_velocity(geometry, radians(d)) / film_h for d in degrees_axis]
    curved_tau = [mu * curved_slide_velocity(geometry, radians(d)) / film_h for d in degrees_axis]

    fig, ax = plt.subplots(figsize=(7.8, 5.2), dpi=140)
    ax.plot(degrees_axis, flat_tau, color="#2e7d32", linewidth=1.8, label="평면(병진)")
    ax.plot(degrees_axis, curved_tau, color="#7a4fb0", linewidth=1.8, label="곡면(진동)")
    ax.axhline(0.0, color="#c7ccd4", linewidth=0.8)
    ax.set_xlim(0, 360)
    ax.set_xticks(range(0, 361, 90))
    ax.set_xlabel("크랭크 각 θ (deg)")
    ax.set_ylabel("전단응력 τ = μU/h (Pa)")
    ax.set_title("면별 전단 트랙션 τ(θ) (§3.6)", fontsize=10, color="none")
    ax.grid(True, color="#e2e5ea", linewidth=0.7)
    _legend_above(ax, ncol=2)
    ax.text(
        0.03,
        0.97,
        f"윤활유 POE VG68, μ = {mu:.3f} Pa·s, 막 10 μm\n"
        f"사이클 평균 마찰 {friction_w * 1.0e3:.0f} mW "
        f"(지시동력 {work.power_w:.0f} W 의 {friction_w / work.power_w * 100.0:.2f} %)\n"
        f"막 누설 {state.total_leakage_m3_s * 1.0e9:.0f} mm³/s",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8.5,
        bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "#c7ccd4"},
    )
    _save(fig, path)


# --------------------------------------------------------------------------
# Figure: gas pressure force and shaft torque
# --------------------------------------------------------------------------


def _moving_average(values: list[float], window: int) -> list[float]:
    """Centred moving average with edge clamping (odd window)."""

    half = window // 2
    count = len(values)
    smoothed = []
    for index in range(count):
        lo = max(0, index - half)
        hi = min(count, index + half + 1)
        smoothed.append(sum(values[lo:hi]) / (hi - lo))
    return smoothed


def _gas_force_sweep(curve: PressureCurve):
    """Return (angles_deg, loads) over the separated range for the gas figures."""

    geometry = curve.geometry
    half_deg = degrees(seal_over_half_angle_rad(geometry)) + 0.25
    angles_deg = [half_deg + (360.0 - 2.0 * half_deg) * i / 400 for i in range(401)]
    loads = [gas_load(geometry, radians(d), trace=curve.trace) for d in angles_deg]
    return angles_deg, loads


def _true_gas_force_sweep(curve: PressureCurve):
    """Return (angles_deg, loads) of the mouth-aware true-geometry gas load.

    Inset an extra half-degree past seal-over so the virtual-work torque's
    central difference never straddles the merge window (where it spikes).
    """

    geometry = curve.geometry
    half_deg = degrees(seal_over_half_angle_rad(geometry)) + 0.5
    angles_deg = [half_deg + (360.0 - 2.0 * half_deg) * i / 400 for i in range(401)]
    loads = [true_gas_load(geometry, radians(d), trace=curve.trace) for d in angles_deg]
    return angles_deg, loads


def render_gas_force_components(curve: PressureCurve, path: Path) -> None:
    """Net gas force components on the rotor over a revolution (PHYSICS.md 4.5)."""

    angles_deg, loads = _gas_force_sweep(curve)
    circ_x = [load.rotor_force_n[0] for load in loads]
    circ_y = [load.rotor_force_n[1] for load in loads]
    true_deg, true_loads = _true_gas_force_sweep(curve)
    true_x = [load.rotor_force_n[0] for load in true_loads]
    true_y = [load.rotor_force_n[1] for load in true_loads]

    fig, ax = plt.subplots(figsize=(8.4, 5.2), dpi=140)
    ax.plot(true_deg, true_x, color=COMPRESSION_COLOR, linewidth=2.0, label="F_x 참 형상")
    ax.plot(true_deg, true_y, color=SUCTION_COLOR, linewidth=2.0, label="F_y 참 형상")
    ax.plot(
        angles_deg, circ_x, color=COMPRESSION_COLOR, linewidth=1.2, linestyle=":", label="F_x 원판"
    )
    ax.plot(angles_deg, circ_y, color=SUCTION_COLOR, linewidth=1.2, linestyle=":", label="F_y 원판")
    ax.axhline(0.0, color="#c7ccd4", linewidth=0.8)
    ax.set_xlabel("크랭크 각 θ (deg)")
    ax.set_ylabel("로터 가스력 F (N)")
    ax.set_title(
        "로터에 작용하는 순 가스력 성분 — 참 형상 vs 원판 (§4.5)", fontsize=10, color="none"
    )
    ax.set_xlim(0, 360)
    ax.set_xticks(range(0, 361, 90))
    ax.grid(True, color="#e2e5ea", linewidth=0.7)
    _legend_above(ax, ncol=2)
    _reed_note(ax)
    _save(fig, path)


def _plot_force_magnitude(curve: PressureCurve, path: Path, *, mode: str) -> None:
    """Net gas force magnitude |F_gas|(theta) on the rotor (= |R_j|, PHYSICS.md 4.5/4.6).

    ``mode`` selects the models drawn: ``"both"`` overlays true (solid) and
    circular (dotted); ``"true"`` and ``"circular"`` draw that single model solid.
    The direction is self-evident (roughly radial through the contact), so the
    magnitude is the reading that matters for the bearing load.
    """

    angles_deg, loads = _gas_force_sweep(curve)
    circ_mag = [hypot(*load.rotor_force_n) for load in loads]
    true_deg, true_loads = _true_gas_force_sweep(curve)
    true_mag = [load.rotor_force_mag_n for load in true_loads]

    fig, ax = plt.subplots(figsize=(8.4, 5.2), dpi=140)
    if mode in ("both", "true"):
        ax.plot(true_deg, true_mag, color=COMPRESSION_COLOR, linewidth=2.0, label="|F_gas| 참 형상")
    if mode in ("both", "circular"):
        ax.plot(
            angles_deg,
            circ_mag,
            color=COMPRESSION_COLOR,
            linewidth=1.2 if mode == "both" else 2.0,
            linestyle=":" if mode == "both" else "-",
            label="|F_gas| 원판",
        )
    peak_deg, peak_mag = (angles_deg, circ_mag) if mode == "circular" else (true_deg, true_mag)
    peak = max(range(len(peak_mag)), key=lambda i: peak_mag[i])
    ax.plot(
        peak_deg[peak],
        peak_mag[peak],
        "o",
        color="#d62728",
        zorder=3,
        label=f"피크 {peak_mag[peak]:.0f} N @ θ={peak_deg[peak]:.0f}°",
    )
    label = {"both": "참 형상 vs 원판", "true": "참 형상", "circular": "원판 근사"}[mode]
    ax.set_xlabel("크랭크 각 θ (deg)")
    ax.set_ylabel("로터 가스력 크기 |F_gas| (N)")
    ax.set_title(f"크랭크 각 vs 가스력 크기 — {label} (§4.5)", fontsize=10, color="none")
    ax.set_xlim(0, 360)
    ax.set_ylim(bottom=0.0)
    ax.set_xticks(range(0, 361, 90))
    ax.grid(True, color="#e2e5ea", linewidth=0.7)
    _legend_above(ax, ncol=3 if mode == "both" else 2)
    _reed_note(ax)
    _save(fig, path)


def render_gas_force_magnitude(curve: PressureCurve, path: Path) -> None:
    """True vs circular gas-force magnitude overlay (PHYSICS.md 4.5)."""

    _plot_force_magnitude(curve, path, mode="both")


def render_gas_force_magnitude_true(curve: PressureCurve, path: Path) -> None:
    """Mouth-aware true-geometry gas-force magnitude alone (PHYSICS.md 4.5)."""

    _plot_force_magnitude(curve, path, mode="true")


def render_gas_force_magnitude_circular(curve: PressureCurve, path: Path) -> None:
    """Circular closed-form gas-force magnitude alone (PHYSICS.md 4.5)."""

    _plot_force_magnitude(curve, path, mode="circular")


def _draw_gas_torque(curve: PressureCurve, path: Path, *, mode: str, write_txt: bool) -> None:
    """Net gas torque about the crank axis, with the pdV cross-check (PHYSICS.md 4.5).

    ``mode`` selects the models drawn: ``"both"`` overlays true (solid) and circular
    (dotted); ``"true"`` and ``"circular"`` draw that single model solid.
    """

    geometry = curve.geometry
    angles_deg, loads = _gas_force_sweep(curve)
    circ_torque = [load.rotor_torque_nm for load in loads]
    true_deg, true_loads = _true_gas_force_sweep(curve)
    # The virtual-work torque differentiates the raster volume trace, so it is
    # stepped by the trace resolution and spikes where the trace's clearance
    # minimum meets seal-over; trim those boundary points and lightly smooth.
    trim = 12
    true_deg = true_deg[trim:-trim]
    true_torque = _moving_average([load.rotor_torque_nm for load in true_loads[trim:-trim]], 9)
    cross = gas_torque_work_j(geometry, trace=curve.trace)
    true_work = true_gas_torque_work_j(geometry, trace=curve.trace)

    fig, ax = plt.subplots(figsize=(8.4, 5.2), dpi=140)
    if mode in ("both", "true"):
        ax.plot(
            true_deg, true_torque, color="#7a4fb0", linewidth=2.0, label="T_gas 참 형상 (가상일)"
        )
    if mode in ("both", "circular"):
        ax.plot(
            angles_deg,
            circ_torque,
            color="#7a4fb0",
            linewidth=1.2 if mode == "both" else 2.0,
            linestyle=":" if mode == "both" else "-",
            label="T_gas = O_r × F (원판)",
        )
    label = {"both": "참 형상 vs 원판", "true": "참 형상 (가상일)", "circular": "원판 (O_r × F)"}[
        mode
    ]
    ax.axhline(0.0, color="#c7ccd4", linewidth=0.8)
    ax.set_xlabel("크랭크 각 θ (deg)")
    ax.set_ylabel("가스 토크 T_gas (N·m)")
    ax.set_title(f"크랭크 축 O 둘레의 순 가스 토크 — {label} (§4.5)", fontsize=10, color="none")
    ax.set_xlim(0, 360)
    ax.set_xticks(range(0, 361, 90))
    ax.grid(True, color="#e2e5ea", linewidth=0.7)
    _legend_above(ax, ncol=1)
    _reed_note(ax)
    _save(fig, path)

    if not write_txt:
        return
    # The former on-figure annotation lives in a sibling .txt (for slides).
    notes = (
        "크랭크 각 - 가스 토크 선도 (gas_force/gas_torque.png)  §4.5\n"
        "========================================================\n\n"
        "축 가스 토크 (= 지시 토크). 베인은 고정체이므로 축 토크에 기여하지\n"
        "않으며, 실링(seal-over) 구간은 적분에서 제외한다.\n\n"
        "두 기준 (같은 압력, 다른 형상):\n"
        "  · 참 형상 (가상일)  T = -(p_c dV_c/dθ + p_s dV_s/dθ)\n"
        "      로터 마우스 캐비티까지 포함한 실제 체적 기준. 부시에 한 일도\n"
        "      포함되므로 O_r x F 로는 재구성되지 않아 가상일로 정의한다.\n"
        "  · 원판 (폐형식)     T = O_r x F\n"
        "      로터 OD 원판 근사. 합력이 O_r 를 지나 모멘트가 닫힌 형태.\n\n"
        "교차검증 (∮ T dθ = -∮ p dV = 지시 일):\n"
        f"  참 형상  ∮ T dθ = {true_work.torque_work_j:.2f} J -> {true_work.power_w:.0f} W"
        f"  (= 지시 동력 {true_work.indicated_power_w:.0f} W, §3.5)\n"
        f"  원판     ∮ T dθ = {cross.torque_work_j:.2f} J -> {cross.power_w:.0f} W\n\n"
        "차이(약 3%, 715->738 W)는 원판이 버리는 로터 마우스 캐비티\n"
        "체적에서 온다. 참 형상 토크가 지시 일과 정확히 일치한다.\n"
    )
    path.with_suffix(".txt").write_text(notes, encoding="utf-8")


def render_gas_torque(curve: PressureCurve, path: Path) -> None:
    """True vs circular gas-torque overlay, with the sibling .txt (PHYSICS.md 4.5)."""

    _draw_gas_torque(curve, path, mode="both", write_txt=True)


def render_gas_torque_true(curve: PressureCurve, path: Path) -> None:
    """Mouth-aware true-geometry gas torque alone (virtual work; PHYSICS.md 4.5)."""

    _draw_gas_torque(curve, path, mode="true", write_txt=False)


def render_gas_torque_circular(curve: PressureCurve, path: Path) -> None:
    """Circular closed-form gas torque alone (O_r × F; PHYSICS.md 4.5)."""

    _draw_gas_torque(curve, path, mode="circular", write_txt=False)


# --------------------------------------------------------------------------
# Figure: crank-pin bearing reaction and drive torque
# --------------------------------------------------------------------------


def _bearing_sweep(curve: PressureCurve):
    """Return (angles_deg, loads) over the separated range for the bearing figures."""

    geometry = curve.geometry
    half_deg = degrees(seal_over_half_angle_rad(geometry)) + 0.25
    angles_deg = [half_deg + (360.0 - 2.0 * half_deg) * i / 400 for i in range(401)]
    loads = [mechanism_load(geometry, radians(d), trace=curve.trace) for d in angles_deg]
    return angles_deg, loads


def _draw_journal_polar(
    curve: PressureCurve, path: Path, *, main: str, overlay_other: bool, write_txt: bool
) -> None:
    """Polar crank-pin bearing-load diagram R_j(theta) = -F_gas (PHYSICS.md 4.6).

    ``main`` selects the coloured locus (``"true"`` mouth-aware or ``"circular"``
    closed form); ``overlay_other`` dots the other model behind it for comparison.
    """

    geometry = curve.geometry
    circ_deg, circ_loads = _bearing_sweep(curve)
    circ_x = [load.journal_force_n[0] for load in circ_loads]
    circ_y = [load.journal_force_n[1] for load in circ_loads]
    circ_mag = [load.journal_force_mag_n for load in circ_loads]

    true_deg, true_loads = _true_gas_force_sweep(curve)
    true_x = [-load.rotor_force_n[0] for load in true_loads]
    true_y = [-load.rotor_force_n[1] for load in true_loads]
    true_mag = [load.rotor_force_mag_n for load in true_loads]

    if main == "true":
        angles_deg, journal_x, journal_y, mags = true_deg, true_x, true_y, true_mag
        other_x, other_y, other_label, model_label = circ_x, circ_y, "원판 폐형식", "참 형상"
    else:
        angles_deg, journal_x, journal_y, mags = circ_deg, circ_x, circ_y, circ_mag
        other_x, other_y, other_label, model_label = true_x, true_y, "참 형상", "원판 근사"
    peak = max(range(len(mags)), key=lambda i: mags[i])

    def _nearest(target_deg: float) -> int:
        # Clamp into the swept range so 0 deg maps to the locus start (~6 deg,
        # where seal-over ends) rather than falling outside it.
        clamped = min(max(target_deg, angles_deg[0]), angles_deg[-1])
        return min(range(len(angles_deg)), key=lambda i: abs(angles_deg[i] - clamped))

    fig, ax = plt.subplots(figsize=(8.2, 6.6), dpi=140)
    if overlay_other:
        ax.plot(
            other_x,
            other_y,
            color="#b8bec8",
            linewidth=1.0,
            linestyle=":",
            zorder=0,
            label=other_label,
        )
    ax.plot(journal_x, journal_y, color="#d7dbe2", linewidth=0.6, zorder=0)
    # Sequential (non-cyclic) map with an absolute 0..360 scale, so each crank
    # angle gets one colour and 0 deg is distinct from 360 deg.
    scatter = ax.scatter(
        journal_x, journal_y, c=angles_deg, cmap="viridis", vmin=0.0, vmax=360.0, s=11, zorder=1
    )
    # Peak reaction vector: the design load, drawn from the crank axis O.
    ax.plot(
        [0.0, journal_x[peak]],
        [0.0, journal_y[peak]],
        color="#d62728",
        linewidth=1.8,
        zorder=2,
        label=f"피크 반력 벡터 ({mags[peak]:.0f} N @ θ={angles_deg[peak]:.0f}°)",
    )
    ax.plot(journal_x[peak], journal_y[peak], "o", color="#d62728", zorder=3)
    ax.plot(0.0, 0.0, "+", color="black", markersize=11, markeredgewidth=1.6, zorder=4)
    ax.text(0.0, 0.0, "  O", fontsize=10, va="center", ha="left", zorder=4)

    # Cardinal reference angles to orient the reader, labelled on the locus.
    for ref_deg in (0.0, 90.0, 180.0, 270.0):
        idx = _nearest(ref_deg)
        ax.plot(
            journal_x[idx],
            journal_y[idx],
            "o",
            markerfacecolor="white",
            markeredgecolor="black",
            markersize=7,
            markeredgewidth=1.3,
            zorder=5,
        )
        label = "θ~0°" if ref_deg == 0.0 else f"θ={ref_deg:.0f}°"
        ax.annotate(
            label,
            xy=(journal_x[idx], journal_y[idx]),
            xytext=(6, 6),
            textcoords="offset points",
            fontsize=9,
            zorder=6,
        )

    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("R_j,x (N, 인렛 +x)")
    ax.set_ylabel("R_j,y (N, 베인 +y)")
    ax.set_title(
        f"저널 베어링 하중 선도  R_j(θ) = -F_gas — {model_label} (§4.6)",
        fontsize=10,
        color="none",
    )
    ax.grid(True, color="#e2e5ea", linewidth=0.7)
    ax.legend(loc="upper left", frameon=False, fontsize=8)
    _reed_note(ax, x=0.02, y=0.02, va="bottom")
    color_bar = fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
    color_bar.set_label("크랭크 각 θ (deg)")
    color_bar.set_ticks(range(0, 361, 45))
    _save(fig, path)

    if not write_txt:
        return
    ports = characteristic_angles(geometry)
    notes = (
        "저널 베어링 하중 선도 (bearing_load/journal_load_polar.png)  §4.6\n"
        "==============================================================\n\n"
        "무엇을 그린 것인가\n"
        "  로터에 작용하는 순 가스력의 반력, 즉 크랭크핀(저널) 베어링이 받는\n"
        "  하중 R_j(θ) = -F_gas(θ) 를 극좌표(벡터 성분)로 나타낸 것.\n"
        "  · 원점 '+O' = 크랭크 축 O.\n"
        "  · 곡선 위의 한 점 = 한 크랭크 각 θ에서의 베어링 하중 벡터. 원점 O에서\n"
        "    그 점까지 그은 화살표가 그 순간의 R_j 이다(방향·크기 모두).\n"
        "  · 가로축 R_j,x = 인렛(+x) 방향, 세로축 R_j,y = 베인(+y) 방향 성분 [N].\n\n"
        "색 (컬러바)\n"
        "  점 색 = 크랭크 각 θ (0->360°, viridis 순차형). 순차형이라 0°(진한 보라)와\n"
        "  360°(노랑)가 서로 다른 색 = 한 바퀴 도는 방향을 색으로 따라갈 수 있다.\n\n"
        "빨간 직선 + 점\n"
        f"  최대 하중이 걸리는 순간의 반력 벡터. 크기 {mags[peak]:.0f} N,\n"
        f"  크랭크 각 θ = {angles_deg[peak]:.0f}° 에서 발생 = 베어링 설계 하중.\n"
        "  직선은 원점 O에서 그 피크 점까지 그은 벡터(방향 = 하중이 향하는 쪽).\n\n"
        "흰 동그라미 (θ = 0 / 90 / 180 / 270°)\n"
        "  방향 판독용 기준각. 곡선을 따라가며 위치를 가늠하는 눈금.\n"
        "  0°는 씰오버 구간이라 실제 점이 없어 곡선 시작점(θ≈6°)에 표시(θ≈0°).\n\n"
        "점선 곡선\n"
        "  원판 폐형식(로터 마우스 무시)의 R_j. 참 형상(색점) 하중이 마우스\n"
        "  캐비티 효과로 그 안쪽 ~20% 로 줄어든다(§4.5).\n\n"
        "참고: 씰오버 구간(|θ| < 약 6°)은 두 챔버가 합쳐져 제외 —\n"
        "곡선이 원점 근처에서 끊긴 이유.\n\n"
        "포트 특성각 (참고, 그림에는 미표시):\n"
        f"  흡입 개시 {degrees(ports.suction_open_rad):.1f}°, "
        f"압축 시작 {degrees(ports.compression_start_rad):.1f}°, "
        f"배출 개시 {degrees(ports.discharge_open_rad):.1f}°, "
        f"배출 종료 {degrees(ports.discharge_close_rad):.1f}°\n"
    )
    path.with_suffix(".txt").write_text(notes, encoding="utf-8")


def render_journal_load_polar(curve: PressureCurve, path: Path) -> None:
    """True bearing-load locus with the circular closed form overlaid (PHYSICS.md 4.6)."""

    _draw_journal_polar(curve, path, main="true", overlay_other=True, write_txt=True)


def render_journal_load_polar_true(curve: PressureCurve, path: Path) -> None:
    """Mouth-aware true-geometry bearing-load locus alone (PHYSICS.md 4.6)."""

    _draw_journal_polar(curve, path, main="true", overlay_other=False, write_txt=False)


def render_journal_load_polar_circular(curve: PressureCurve, path: Path) -> None:
    """Circular closed-form bearing-load locus alone (PHYSICS.md 4.6)."""

    _draw_journal_polar(curve, path, main="circular", overlay_other=False, write_txt=False)


def _draw_drive_torque(curve: PressureCurve, path: Path, *, basis: str) -> None:
    """Drive torque vs indicated torque with the shaft-power balance (PHYSICS.md 4.6/4.7).

    ``basis`` picks the plotted gas torque: ``"circular"`` (closed form, 715 W) or
    ``"true"`` (mouth-aware virtual-work torque, 738 W). Friction is the same for
    both (velocity-driven bush + load-independent Petroff journal), so only the
    ``T_gas``/``T_drive`` curves shift between the two.
    """

    geometry = curve.geometry
    work = shaft_work_j(geometry, trace=curve.trace)
    if basis == "circular":
        angles_deg, loads = _bearing_sweep(curve)
        gas_torque = [load.gas_torque_nm for load in loads]
        drive_torque = [load.drive_torque_nm for load in loads]
        tgas_power = gas_torque_work_j(geometry, trace=curve.trace).power_w
        basis_label = "원판"
    else:
        true_deg, true_loads = _true_gas_force_sweep(curve)
        friction = [
            mechanism_load(geometry, radians(d), trace=curve.trace).friction_torque_nm
            for d in true_deg
        ]
        gas_raw = [load.rotor_torque_nm for load in true_loads]
        # The virtual-work torque spikes at the seal-over ends and is stepped by
        # the trace resolution; trim those points and lightly smooth (as in the
        # gas-torque figure), then add the smooth friction.
        trim = 12
        angles_deg = true_deg[trim:-trim]
        gas_torque = _moving_average(gas_raw[trim:-trim], 9)
        drive_torque = [g + f for g, f in zip(gas_torque, friction[trim:-trim], strict=True)]
        tgas_power = true_gas_torque_work_j(geometry, trace=curve.trace).power_w
        basis_label = "참 형상"

    fig, ax = plt.subplots(figsize=(8.4, 5.4), dpi=140)
    ax.plot(
        angles_deg,
        gas_torque,
        color=COMPRESSION_COLOR,
        linewidth=2.0,
        label=f"T_gas (지시, {basis_label})",
    )
    ax.plot(
        angles_deg,
        drive_torque,
        color="#7a4fb0",
        linewidth=1.4,
        linestyle="--",
        label="T_drive = T_gas + 부시+저널 마찰",
    )
    ax.axhline(0.0, color="#c7ccd4", linewidth=0.8)
    ax.set_xlabel("크랭크 각 θ (deg)")
    ax.set_ylabel("토크 (N·m)")
    ax.set_title(f"구동 토크와 지시 토크 — {basis_label} (§4.6·§4.7)", fontsize=10, color="none")
    ax.set_xlim(0, 360)
    ax.set_xticks(range(0, 361, 90))
    ax.grid(True, color="#e2e5ea", linewidth=0.7)
    _legend_above(ax, ncol=2)
    ax.text(
        0.03,
        0.55,
        f"축동력 {work.shaft_power_w:.0f} W = 이상지시 {work.indicated_power_w:.0f}"
        f" + 부시 {work.bush_friction_power_w:.2f} + 저널 {work.journal_friction_power_w:.1f} W\n"
        "리드밸브 미포함 (§3.8 과압은 별도 성능 항 · loads 미반영)\n"
        f"피크 베어링 하중 {work.peak_journal_force_n:.0f} N\n"
        f"표시 T_gas: {basis_label} 기준 (∮T·f = {tgas_power:.0f} W)",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "#c7ccd4"},
    )
    _save(fig, path)


def render_drive_torque(curve: PressureCurve, path: Path) -> None:
    """Drive torque on the mouth-aware true gas torque (headline, PHYSICS.md 4.6/4.7)."""

    _draw_drive_torque(curve, path, basis="true")


def render_drive_torque_circular(curve: PressureCurve, path: Path) -> None:
    """Drive torque on the circular closed-form gas torque (PHYSICS.md 4.6/4.7)."""

    _draw_drive_torque(curve, path, basis="circular")


def _draw_friction_breakdown(curve: PressureCurve, path: Path, *, basis_label: str) -> None:
    """Journal vs bush-film friction power over a revolution (PHYSICS.md 3.6/4.7).

    The two friction terms are **load-independent** — velocity-driven bush film and
    Petroff (concentric) journal — so the curves are identical whether the gas
    force is the circular closed form or the mouth-aware true geometry;
    ``basis_label`` only tags the (hidden) title. A small note on the axes states
    this so the paired true/circular slides read honestly.
    """

    geometry = curve.geometry
    shaft_speed = geometry.angular_speed_rad_s
    angles_deg, loads = _bearing_sweep(curve)
    bush_power = [load.bush_friction_torque_nm * shaft_speed for load in loads]
    journal_power = [load.journal_friction_torque_nm * shaft_speed for load in loads]
    work = shaft_work_j(geometry, trace=curve.trace)

    fig, ax = plt.subplots(figsize=(8.4, 5.4), dpi=140)
    ax.plot(
        angles_deg,
        journal_power,
        color="#c0392b",
        linewidth=2.0,
        label=f"저널 (Petroff, 평균 {work.journal_friction_power_w:.1f} W)",
    )
    ax.plot(
        angles_deg,
        bush_power,
        color="#2e7d32",
        linewidth=2.0,
        label=f"부시막 (평균 {work.bush_friction_power_w:.2f} W)",
    )
    ax.axhline(0.0, color="#c7ccd4", linewidth=0.8)
    ax.set_xlabel("크랭크 각 θ (deg)")
    ax.set_ylabel("마찰 손실 (W)")
    ax.set_title(
        f"마찰 손실 분해 — 저널이 지배 · {basis_label} (§3.6·§4.7)", fontsize=10, color="none"
    )
    ax.set_xlim(0, 360)
    ax.set_xticks(range(0, 361, 90))
    ax.grid(True, color="#e2e5ea", linewidth=0.7)
    ax.text(
        0.03,
        0.97,
        "마찰은 하중 독립적 (속도구동 부시 + Petroff 저널)\n"
        "→ 원판·참 형상 basis 무관, 두 버전 동일",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8,
        color="#555",
        bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "#d7dbe2"},
    )
    _reed_note(ax, x=0.98, y=0.02, ha="right", va="bottom")
    _legend_above(ax, ncol=2)
    _save(fig, path)


def render_friction_breakdown_true(curve: PressureCurve, path: Path) -> None:
    """Friction breakdown, true-geometry slide (identical curves; PHYSICS.md 3.6/4.7)."""

    _draw_friction_breakdown(curve, path, basis_label="참 형상")


def render_friction_breakdown_circular(curve: PressureCurve, path: Path) -> None:
    """Friction breakdown, circular slide (identical curves; PHYSICS.md 3.6/4.7)."""

    _draw_friction_breakdown(curve, path, basis_label="원판")


def render_bush_uniform_clearance_model(geometry: RotaryGeometry, path: Path) -> None:
    """Schematic of the swing-bush uniform-clearance film friction model (PHYSICS.md 3.6).

    Left: where the two 10 um films sit (bush flat vs vane, bush curve vs groove).
    Right: the Couette idealisation — a uniform gap h sheared at the sliding speed
    U, tau = mu U / h. Pedagogical; clearances are exaggerated (not to scale).
    """

    bush = SwingBush()
    flat_h, curved_h = film_thicknesses_m(geometry, bush)
    mu = LUBRICANT_VISCOSITY_PA_S
    state = prescribed_state(geometry, 0.0)
    gx_m, gy_m = state.cutout_center_m
    gx, gy = gx_m / MM, gy_m / MM
    cut_r = geometry.cutout_radius_m / MM
    half_w = 0.5 * geometry.vane_width_m / MM

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13.6, 6.2), dpi=140)
    fig.suptitle("스윙 부시 균일 클리어런스 유막 마찰 모델 (§3.6)", fontsize=13, color="none")

    # ---- Panel A: real placement (zoom on the groove) ----
    ax_a.set_aspect("equal")
    ax_a.axis("off")
    ax_a.set_xlim(gx - 13, gx + 13)
    ax_a.set_ylim(gy - 11, gy + 13.5)
    ax_a.add_patch(
        plt.Circle((gx, gy), cut_r, facecolor="#eef1f5", edgecolor="#888", linewidth=1.4, zorder=1)
    )
    ax_a.add_patch(
        Rectangle(
            (-half_w, gy - 10.5),
            2 * half_w,
            23.0,
            facecolor=VANE_COLOR,
            edgecolor="black",
            linewidth=1.0,
            zorder=2,
        )
    )
    for side in (1.0, -1.0):
        xs, ys = _bush_outline_mm(gx_m + side * bush.piece_shift_m, gy_m, side)
        ax_a.fill(xs, ys, facecolor="#9fb6d4", edgecolor="black", linewidth=1.0, zorder=3)
    ax_a.annotate(
        f"평면 유막 h ~ {flat_h * 1e6:.0f} " + r"$\mu$m" + "\n(부시 평면 ↔ 베인)",
        xy=(half_w, gy - 2.5),
        xytext=(gx + 4.5, gy - 9.5),
        fontsize=9.5,
        color="#c0392b",
        arrowprops={"arrowstyle": "->", "color": "#c0392b"},
    )
    ax_a.annotate(
        f"곡면 유막 h ~ {curved_h * 1e6:.0f} " + r"$\mu$m" + "\n(부시 곡면 ↔ 그루브)",
        xy=(gx - cut_r + 0.3, gy + 2.0),
        xytext=(gx - 12.5, gy + 9.0),
        fontsize=9.5,
        color="#2e7d32",
        arrowprops={"arrowstyle": "->", "color": "#2e7d32"},
    )
    ax_a.text(
        gx, gy + 12.4, "실제 배치 (유막 10 " + r"$\mu$m — 육안 미세)", ha="center", fontsize=10
    )

    # ---- Panel B: Couette idealisation ----
    ax_b.set_xlim(0, 10)
    ax_b.set_ylim(0, 10)
    ax_b.axis("off")
    x_l, x_r, y_b, y_t = 1.6, 8.4, 3.3, 5.9
    ax_b.add_patch(Rectangle((x_l, y_b), x_r - x_l, y_t - y_b, facecolor="#fde8c8", zorder=0))
    ax_b.add_patch(
        Rectangle((x_l, y_b - 0.6), x_r - x_l, 0.6, facecolor="#9fb6d4", edgecolor="black")
    )
    ax_b.add_patch(Rectangle((x_l, y_t), x_r - x_l, 0.6, facecolor=VANE_COLOR, edgecolor="black"))
    for i in range(7):
        frac = i / 6.0
        yy = y_b + (y_t - y_b) * frac
        ax_b.annotate(
            "",
            xy=(x_l + 1.0 + 3.1 * frac, yy),
            xytext=(x_l + 1.0, yy),
            arrowprops={"arrowstyle": "->", "color": "#7a4fb0", "lw": 1.2},
        )
    ax_b.annotate(
        "",
        xy=(x_r - 0.6, y_t + 0.3),
        xytext=(x_r - 2.6, y_t + 0.3),
        arrowprops={"arrowstyle": "-|>", "color": "black", "lw": 2.0},
    )
    ax_b.text(x_r - 0.4, y_t + 0.3, "U", va="center", fontsize=13)
    ax_b.annotate(
        "",
        xy=(x_l - 0.45, y_t),
        xytext=(x_l - 0.45, y_b),
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )
    ax_b.text(x_l - 0.75, 0.5 * (y_b + y_t), "h", ha="right", va="center", fontsize=13)
    ax_b.text(
        0.5 * (x_l + x_r),
        y_t + 1.2,
        "이동벽 (베인·그루브 상대 미끄럼 U)",
        ha="center",
        fontsize=9.5,
    )
    ax_b.text(
        0.5 * (x_l + x_r), y_b - 0.3, "부시 표면 (고정)", ha="center", va="center", fontsize=9
    )
    ax_b.text(5, 9.2, "모델: 균일 유막 Couette 전단", ha="center", fontsize=11)
    ax_b.text(5, 1.9, r"$\tau = \mu\,U/h$", ha="center", fontsize=15)
    ax_b.text(
        5,
        0.7,
        rf"$\mu$={mu:.3f} Pa·s (POE VG68),  h ~ 10 "
        + r"$\mu$m 균일"
        + "\nU ≤ 0.86 m/s (평면), ≤ 0.27 (곡면)",
        ha="center",
        fontsize=9,
        color="#555",
    )
    _save(fig, path)


def render_journal_concentric_clearance_model(geometry: RotaryGeometry, path: Path) -> None:
    """Schematic of the crank-pin journal concentric-clearance (Petroff) model (PHYSICS.md 4.7).

    Left: cross-section — the crank pin sits **concentrically** in the rotor bore
    with a uniform radial clearance c_j (the Petroff assumption, eps=0). Right:
    axial section + Petroff law. Clearances are exaggerated (real c_j = 15 um).
    """

    r_j = JOURNAL_RADIUS_M / MM
    l_j = JOURNAL_LENGTH_M / MM
    c_j_um = JOURNAL_CLEARANCE_M * 1e6
    rotor_r = geometry.rotor_radius_m / MM
    c_draw = 0.09 * r_j  # exaggerated radial clearance for drawing

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13.6, 6.4), dpi=140)
    fig.suptitle("크랭크핀 저널 동심 클리어런스 마찰 모델 (§4.7)", fontsize=13, color="none")

    # ---- Panel A: cross-section (concentric pin in the rotor bore) ----
    ax_a.set_aspect("equal")
    ax_a.axis("off")
    lim = rotor_r + 3.0
    ax_a.set_xlim(-lim, lim)
    ax_a.set_ylim(-lim, lim + 2.0)
    ax_a.add_patch(
        plt.Circle(
            (0, 0), rotor_r, facecolor=ROTOR_COLOR, edgecolor="black", linewidth=1.8, zorder=1
        )
    )
    ax_a.add_patch(
        plt.Circle(
            (0, 0), r_j + c_draw, facecolor="#fde8c8", edgecolor="#c0392b", linewidth=1.2, zorder=2
        )
    )
    ax_a.add_patch(
        plt.Circle((0, 0), r_j, facecolor="#6b7a8f", edgecolor="black", linewidth=1.5, zorder=3)
    )
    ax_a.plot(0, 0, "+", color="white", markersize=9, markeredgewidth=1.6, zorder=4)
    ax_a.annotate(
        "",
        xy=(0.6 * r_j * cos(radians(30)), 0.6 * r_j * sin(radians(30))),
        xytext=(0.6 * r_j * cos(radians(150)), 0.6 * r_j * sin(radians(150))),
        arrowprops={
            "arrowstyle": "-|>",
            "color": "white",
            "lw": 2.0,
            "connectionstyle": "arc3,rad=-0.4",
        },
        zorder=5,
    )
    ax_a.text(
        0,
        -0.32 * r_j,
        "크랭크핀\n" + r"$\omega_{rel}$",
        ha="center",
        va="center",
        color="white",
        fontsize=11,
        zorder=6,
    )
    ax_a.annotate(
        "",
        xy=(r_j * cos(radians(-22)), r_j * sin(radians(-22))),
        xytext=(0, 0),
        arrowprops={"arrowstyle": "->", "color": "black"},
        zorder=6,
    )
    ax_a.text(
        0.5 * r_j * cos(radians(-22)) + 0.6,
        0.5 * r_j * sin(radians(-22)) - 1.4,
        r"$r_j$=14.2",
        fontsize=11,
        zorder=6,
    )
    ax_a.annotate(
        r"$c_j$ 균일 (동심)" + f"\n실제 {c_j_um:.0f} " + r"$\mu$m 과장",
        xy=(0, r_j + 0.5 * c_draw),
        xytext=(-0.55 * lim, 0.6 * rotor_r),
        fontsize=10,
        color="#c0392b",
        arrowprops={"arrowstyle": "->", "color": "#c0392b"},
    )
    ax_a.text(
        0,
        rotor_r + 1.4,
        "단면 — 동심 클리어런스 (Petroff, " + r"$\varepsilon$=0)",
        ha="center",
        fontsize=11,
    )
    ax_a.text(0.6 * rotor_r, -0.66 * rotor_r, "로터", fontsize=11)

    # ---- Panel B: axial section + Petroff law ----
    ax_b.set_xlim(0, 10)
    ax_b.set_ylim(0, 10)
    ax_b.axis("off")
    x0, x1, yc, ph, gap = 1.8, 8.0, 6.7, 2.0, 0.45
    ax_b.add_patch(Rectangle((x0, yc + ph), x1 - x0, gap, facecolor="#fde8c8", zorder=0))
    ax_b.add_patch(Rectangle((x0, yc - ph - gap), x1 - x0, gap, facecolor="#fde8c8", zorder=0))
    ax_b.add_patch(
        Rectangle((x0, yc - ph), x1 - x0, 2 * ph, facecolor="#6b7a8f", edgecolor="black")
    )
    ax_b.plot([x0, x1], [yc + ph + gap, yc + ph + gap], color="#c0392b", lw=1.2)
    ax_b.plot([x0, x1], [yc - ph - gap, yc - ph - gap], color="#c0392b", lw=1.2)
    ax_b.annotate(
        "",
        xy=(x1, yc - ph - gap - 0.6),
        xytext=(x0, yc - ph - gap - 0.6),
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )
    ax_b.text(
        0.5 * (x0 + x1), yc - ph - gap - 1.3, rf"$L_j$ = {l_j:.0f} mm", ha="center", fontsize=11
    )
    ax_b.annotate(
        "",
        xy=(x1 + 0.55, yc + ph),
        xytext=(x1 + 0.55, yc - ph),
        arrowprops={"arrowstyle": "<->", "color": "black"},
    )
    ax_b.text(x1 + 0.85, yc, r"$2r_j$", va="center", fontsize=11)
    ax_b.annotate(
        r"$c_j$ (유막, 균일)",
        xy=(x1 - 1.2, yc + ph + 0.5 * gap),
        xytext=(x1 - 2.6, yc + ph + 1.5),
        fontsize=10,
        color="#c0392b",
        arrowprops={"arrowstyle": "->", "color": "#c0392b"},
    )
    ax_b.text(2.6, 9.5, "축방향 단면", ha="center", fontsize=11)
    ax_b.text(
        5,
        1.9,
        r"$T_j = \dfrac{2\pi\,\mu\,|\omega_{rel}|\,r_j^{3}\,L_j}{c_j}$",
        ha="center",
        fontsize=15,
    )
    ax_b.text(
        5,
        0.5,
        r"동심·하중 독립 (Petroff).  $r_j/c_j$~947,  $\mu$=0.010 Pa·s"
        + "\n실제는 편심 (S~0.07) → 마찰 과소평가",
        ha="center",
        fontsize=9,
        color="#555",
    )
    _save(fig, path)


def render_journal_relative_speed(geometry: RotaryGeometry, path: Path) -> None:
    """Crank-pin <-> rotor relative angular speed omega_rel over a revolution (PHYSICS.md 4.7).

    The crank pin turns at the constant shaft speed while the rotor only swings
    (+/-10.37 deg), so its spin oscillates; their difference is the oil-shear speed
    the Petroff journal friction uses, ``omega_rel = omega (1 - dphi/dtheta)``.
    """

    omega = geometry.angular_speed_rad_s
    n = 361
    angles_deg = [360.0 * i / (n - 1) for i in range(n)]
    omega_rel = [journal_relative_speed_rad_s(geometry, radians(d)) for d in angles_deg]
    omega_rotor = [omega - w for w in omega_rel]
    lo, hi = min(omega_rel), max(omega_rel)

    fig, ax = plt.subplots(figsize=(8.4, 5.2), dpi=140)
    ax.plot(
        angles_deg,
        omega_rel,
        color="#c0392b",
        linewidth=2.2,
        label=r"$\omega_{rel}=\omega(1-\mathrm{d}\phi/\mathrm{d}\theta)$ (크랭크핀↔로터)",
    )
    ax.axhline(
        omega,
        color="#7a4fb0",
        linestyle="--",
        linewidth=1.4,
        label=rf"$\omega$ 샤프트 (일정 {omega:.0f} rad/s)",
    )
    ax.plot(
        angles_deg,
        omega_rotor,
        color="#2e7d32",
        linewidth=1.4,
        label=r"$\omega_{rotor}=\omega\,\mathrm{d}\phi/\mathrm{d}\theta$  (로터 스핀, ±스윙)",
    )
    ax.axhline(0.0, color="#c7ccd4", linewidth=0.8)
    ax.set_xlabel("크랭크 각 θ (deg)")
    ax.set_ylabel("각속도 (rad/s)")
    ax.set_title("샤프트–로터 상대 각속도 — 저널 오일 전단 속도 (§4.7)", fontsize=10, color="none")
    ax.set_xlim(0, 360)
    ax.set_xticks(range(0, 361, 90))
    ax.grid(True, color="#e2e5ea", linewidth=0.7)
    ax.text(
        0.03,
        0.60,
        rf"$\omega_{{rel}}$ ~ {lo:.0f}–{hi:.0f} rad/s (샤프트 {omega:.0f} $\pm$ 스윙)"
        + "\nPetroff 마찰 "
        + r"$T_j \propto |\omega_{rel}|$"
        + " → 저널 손실 8.9 W (§4.7)",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "#c7ccd4"},
    )
    _legend_above(ax, ncol=1)
    _save(fig, path)


def render_journal_eccentricity(curve: PressureCurve, path: Path) -> None:
    """Crank-pin journal eccentricity ratio and minimum oil film over a revolution (PHYSICS.md 4.9).

    The short-bearing (Ocvirk) film is balanced against the true gas reaction
    ``|R_j(theta)|`` (Section 4.6) at the crank-pin entrainment speed, giving the
    running eccentricity ``eps = e/c_j`` and the minimum film ``c_j(1-eps)``. The
    peak eccentricity ``eps ~ 0.71`` (minimum film ~4.4 um from the 15 um
    clearance) sits at the peak reaction near 226 deg -- the concentric Petroff
    model of Section 4.7 ignores this.
    """

    cycle = eccentricity_cycle(curve.geometry, trace=curve.trace)
    theta = [degrees(a) for a in cycle.crank_angle_rad]
    eps = list(cycle.eccentricity_ratio)
    film_um = [h * 1.0e6 for h in cycle.min_film_thickness_m]
    ipeak = max(range(len(eps)), key=lambda k: eps[k])

    fig, ax = plt.subplots(figsize=(8.4, 5.2), dpi=140)
    (line_eps,) = ax.plot(
        theta, eps, color="#c0392b", linewidth=2.3, label=r"편심비 $\epsilon=e/c_j$"
    )
    ax.axhline(1.0, color="#c7ccd4", linestyle=":", linewidth=1.0)
    ax.set_xlabel("크랭크 각 θ (deg)")
    ax.set_ylabel(r"편심비 $\epsilon$", color="#c0392b")
    ax.set_ylim(0.0, 1.0)
    ax.set_xlim(0, 360)
    ax.set_xticks(range(0, 361, 90))
    ax.tick_params(axis="y", labelcolor="#c0392b")
    ax.grid(True, color="#e2e5ea", linewidth=0.7)
    ax.set_title("크랭크핀 저널 편심비·최소 유막 — 단축 베어링 (§4.9)", fontsize=10, color="none")

    ax2 = ax.twinx()
    (line_film,) = ax2.plot(
        theta,
        film_um,
        color="#2e7d32",
        linewidth=1.9,
        label=r"최소 유막 $h_{min}=c_j(1-\epsilon)$",
    )
    ax2.set_ylabel(r"최소 유막 $h_{min}$ (µm)", color="#2e7d32")
    ax2.tick_params(axis="y", labelcolor="#2e7d32")
    ax2.set_ylim(0.0, JOURNAL_CLEARANCE_M * 1.0e6)

    ax.annotate(
        rf"peak $\epsilon$={eps[ipeak]:.2f} @ {theta[ipeak]:.0f}°"
        + f"\n$h_{{min}}$={film_um[ipeak]:.1f} µm (|$R_j$|={cycle.load_n[ipeak] / 1e3:.2f} kN)",
        xy=(theta[ipeak], eps[ipeak]),
        xytext=(theta[ipeak] - 150.0, 0.86),
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "#c7ccd4"},
        arrowprops={"arrowstyle": "->", "color": "#8a8f98"},
    )
    ax.text(
        0.03,
        0.06,
        rf"$\epsilon$ {min(eps):.2f}–{max(eps):.2f}, $c_j$={JOURNAL_CLEARANCE_M * 1e6:.0f} µm 가정"
        + "\n하중 의존 유막 (Petroff 동심 모델 대비)",
        transform=ax.transAxes,
        va="bottom",
        ha="left",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "#c7ccd4"},
    )
    ax.legend(
        handles=[line_eps, line_film],
        loc="lower center",
        bbox_to_anchor=(0.5, 1.0),
        ncol=2,
        fontsize=9,
        framealpha=0.9,
    )
    _save(fig, path)


def render_eccentric_friction_power(curve: PressureCurve, path: Path) -> None:
    """Eccentric vs concentric (Petroff) crank-pin friction over a revolution (PHYSICS.md 4.9).

    The short-bearing Couette friction ``T_petroff / sqrt(1-eps^2)`` grows with the
    load-driven eccentricity, so it rises above the load-independent Petroff
    estimate (Section 4.7) where the reaction peaks -- lifting the cycle-mean
    journal loss from ~9.2 W to ~10.8 W (~18%).
    """

    cycle = eccentricity_cycle(curve.geometry, trace=curve.trace)
    theta = [degrees(a) for a in cycle.crank_angle_rad]

    fig, ax = plt.subplots(figsize=(8.4, 5.2), dpi=140)
    ax.plot(
        theta,
        cycle.friction_power_w,
        color="#c0392b",
        linewidth=2.3,
        label=r"편심 마찰 $P_f=T_f|\omega_{rel}|,\; T_f=T_{Petroff}/\sqrt{1-\epsilon^2}$",
    )
    # Concentric Petroff power at each angle: P_f * sqrt(1-eps^2) (removes the
    # eccentricity factor, recovering the load-independent shear of Section 4.7).
    petroff_curve = [
        p * (1.0 - e**2) ** 0.5
        for p, e in zip(cycle.friction_power_w, cycle.eccentricity_ratio, strict=True)
    ]
    ax.plot(
        theta,
        petroff_curve,
        color="#2e7d32",
        linewidth=1.7,
        linestyle="--",
        label=r"Petroff 동심 $\epsilon=0$ (§4.7)",
    )
    ax.axhline(
        cycle.mean_friction_power_w,
        color="#c0392b",
        linestyle=":",
        linewidth=1.2,
        label=f"편심 평균 {cycle.mean_friction_power_w:.1f} W",
    )
    ax.axhline(
        cycle.petroff_mean_friction_power_w,
        color="#2e7d32",
        linestyle=":",
        linewidth=1.2,
        label=f"Petroff 평균 {cycle.petroff_mean_friction_power_w:.1f} W",
    )
    ax.set_xlabel("크랭크 각 θ (deg)")
    ax.set_ylabel("저널 마찰 손실 (W)")
    ax.set_xlim(0, 360)
    ax.set_xticks(range(0, 361, 90))
    ax.set_ylim(bottom=0.0)
    ax.grid(True, color="#e2e5ea", linewidth=0.7)
    ax.set_title("편심 저널 마찰 vs Petroff 동심 (§4.9)", fontsize=10, color="none")
    _legend_above(ax, ncol=2)
    _save(fig, path)


def render_journal_film_coordinates(geometry: RotaryGeometry, path: Path) -> None:
    """Rotor eccentricity and crank-pin journal film coordinates (PHYSICS.md 4.10).

    Panel A: the full rotor, eccentric from the shaft axis O by e=4.5 mm (crank
    throw); the crank-pin bore sits at the rotor centre O_r=O_b. Panel B: the
    journal zoom (clearance ~300x exaggerated) with the 2D bearing eccentricity
    e = O_b - O_j (<= c_j ~ 15 um, the EOM state), angle alpha, and film h(alpha).
    """

    theta = radians(45.0)
    state = prescribed_state(geometry, theta)
    rx, ry = state.rotor_center_m[0] / MM, state.rotor_center_m[1] / MM
    bore = geometry.cylinder_radius_m / MM
    r_j = JOURNAL_RADIUS_M / MM

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14.0, 7.0), dpi=140)
    fig.suptitle("로터 편심과 크랭크핀 저널 좌표 (§4.10)", fontsize=13, color="none")

    # ---- Panel A: the full rotor, macro 4.5 mm eccentricity ----
    ax_a.set_aspect("equal")
    ax_a.axis("off")
    lim = bore + 6
    ax_a.set_xlim(-lim, lim)
    ax_a.set_ylim(-lim, lim)
    ax_a.add_patch(
        plt.Circle((0, 0), bore, facecolor="#f4f6f9", edgecolor="black", lw=2.0, zorder=1)
    )
    contour = rotor_contour(geometry, theta)
    ax_a.fill(
        [p[0] / MM for p in contour.material],
        [p[1] / MM for p in contour.material],
        facecolor=ROTOR_COLOR,
        edgecolor="black",
        lw=1.5,
        zorder=2,
    )
    bush = SwingBush()
    gx, gy = state.cutout_center_m
    for side in (1.0, -1.0):
        xs, ys = _bush_outline_mm(gx + side * bush.piece_shift_m, gy, side)
        ax_a.fill(xs, ys, facecolor=BUSH_COLOR, edgecolor="black", lw=0.8, zorder=3)
    ax_a.add_patch(
        plt.Circle((rx, ry), r_j, facecolor="#dbe7f5", edgecolor="#2166ac", lw=1.6, zorder=4)
    )
    ax_a.add_patch(
        plt.Circle((rx, ry), 3.0, facecolor="#6b7a8f", edgecolor="black", lw=1.0, zorder=5)
    )
    ax_a.plot(0, 0, "+", color="black", ms=14, mew=2.2, zorder=6)
    ax_a.text(-1.5, -3.0, "O (샤프트축 = 실린더 중심)", fontsize=10, ha="right", zorder=6)
    ax_a.annotate(
        "",
        xy=(rx, ry),
        xytext=(0, 0),
        arrowprops={"arrowstyle": "-|>", "color": "#c0392b", "lw": 2.2},
        zorder=6,
    )
    ax_a.text(
        rx * 0.5 - 11, ry * 0.5 + 1.5, "e = 4.5 mm\n(크랭크 throw)", color="#c0392b", fontsize=10
    )
    ax_a.plot(rx, ry, "o", color="#c0392b", ms=6, zorder=7)
    ax_a.text(
        rx + 1.2,
        ry + 1.4,
        r"$O_r=O_b$ (로터 중심)" + "\n" + r"$\approx O_j$ (크랭크핀)",
        fontsize=9.5,
        zorder=7,
    )
    ax_a.add_patch(
        plt.Circle((rx, ry), r_j + 1.5, fill=False, ec="#c0392b", lw=1.2, ls="--", zorder=6)
    )
    ax_a.annotate(
        "확대 → B",
        xy=(rx + r_j, ry),
        xytext=(rx + r_j + 2, ry - 6),
        fontsize=10,
        color="#c0392b",
        arrowprops={"arrowstyle": "->", "color": "#c0392b"},
    )
    ax_a.text(
        -lim + 1, lim - 3, "A. 전체 로터 — 매크로 편심 e=4.5mm", fontsize=12, fontweight="bold"
    )

    # ---- Panel B: journal zoom, e = O_b - O_j (clearance exaggerated) ----
    r_bz, r_jz = 1.0, 0.72
    ehat = 0.6 * (r_bz - r_jz)
    ang = radians(52)
    ex, ey = ehat * cos(ang), ehat * sin(ang)  # e = O_b - O_j, O_b at origin
    ojx, ojy = -ex, -ey  # O_j = O_b - e
    ring = [(cos(radians(t)), sin(radians(t))) for t in range(0, 361, 3)]
    ax_b.fill([r_bz * x for x, _ in ring], [r_bz * y for _, y in ring], color="#eef4fb", zorder=0)
    ax_b.plot(
        [r_bz * x for x, _ in ring], [r_bz * y for _, y in ring], "--", color="#8a8f98", lw=1.6
    )
    ax_b.add_patch(plt.Circle((ojx, ojy), r_jz, fill=False, ec="#2166ac", lw=2.3, zorder=2))

    ax_b.plot(0, 0, "o", color="#c0392b", ms=6, zorder=4)
    ax_b.annotate(r"$O_b$ (로터 보어=슬리브)", (0, 0), (0.06, 0.12), fontsize=10, color="#c0392b")
    ax_b.plot(ojx, ojy, "o", color="#2166ac", ms=6, zorder=4)
    ax_b.annotate(
        r"$O_j$ (크랭크핀=저널)", (ojx, ojy), (ojx - 0.1, ojy - 0.22), fontsize=10, color="#2166ac"
    )

    # eccentricity vector e = O_b - O_j (from O_j to O_b)
    ax_b.annotate(
        "",
        xy=(0, 0),
        xytext=(ojx, ojy),
        arrowprops={"arrowstyle": "-|>", "color": "#c0392b", "lw": 2.0},
    )
    ax_b.annotate(
        r"$\vec e=O_b-O_j$",
        (ojx / 2, ojy / 2),
        (ojx / 2 + 0.1, ojy / 2 - 0.02),
        color="#c0392b",
        fontsize=11,
    )

    # global-parallel axes at O_b
    for xy, txt, off in (
        ((0.6, 0), r"$\hat x$", (0.5, -0.13)),
        ((0, 0.6), r"$\hat y$", (0.04, 0.54)),
    ):
        ax_b.annotate("", xy=xy, xytext=(0, 0), arrowprops={"arrowstyle": "-|>", "color": "black"})
        ax_b.annotate(txt, xy, off, fontsize=10)

    # circumferential angle alpha + film thickness h(alpha) (thick side)
    a_s = radians(60)
    ax_b.plot([0, r_bz * cos(a_s)], [0, r_bz * sin(a_s)], color="#2e7d32", lw=1.2)
    ax_b.add_patch(Arc((0, 0), 0.66, 0.66, angle=0, theta1=0, theta2=60, color="#2e7d32", lw=1.5))
    ax_b.annotate(r"$\alpha$", (0.36, 0.22), (0.36, 0.22), color="#2e7d32", fontsize=12)
    bx, by = r_bz * cos(a_s), r_bz * sin(a_s)
    jx, jy = ojx + r_jz * cos(a_s), ojy + r_jz * sin(a_s)
    ax_b.annotate(
        "",
        xy=(bx, by),
        xytext=(jx, jy),
        arrowprops={"arrowstyle": "<|-|>", "color": "#7a4fb0", "lw": 1.6},
    )
    ax_b.annotate(r"$h(\alpha)$", (bx, by), (bx + 0.02, by + 0.05), color="#7a4fb0", fontsize=11)
    mx, my = -cos(ang), -sin(ang)
    ax_b.annotate(
        r"$h_{min}$",
        (r_bz * mx, r_bz * my),
        (r_bz * mx - 0.05, r_bz * my - 0.18),
        color="#c0392b",
        fontsize=10,
    )

    ax_b.text(
        -1.32,
        1.34,
        "B. 저널 확대 — 2D 편심 (클리어런스 300× 과장)",
        fontsize=11.5,
        fontweight="bold",
    )
    ax_b.text(
        -1.34,
        -1.52,
        r"$\vec e=O_b-O_j=(e_x,e_y),\ \ |\vec e|\leq c_j\approx15\,\mu$m"
        "\n"
        r"$h(\alpha)=c_j+(e_x\cos\alpha+e_y\sin\alpha)$",
        fontsize=11,
        bbox={"boxstyle": "round", "fc": "white", "ec": "#c7ccd4"},
    )
    ax_b.set_aspect("equal")
    ax_b.axis("off")
    ax_b.set_xlim(-1.4, 1.4)
    ax_b.set_ylim(-1.65, 1.5)

    _save(fig, path)


def render_journal_film_axial(geometry: RotaryGeometry, path: Path) -> None:
    """Journal film axial section — tilt / 2D Reynolds extension (PHYSICS.md 4.10).

    Axial coordinate zeta in [-L_j/2, L_j/2]; journal-axis tilt tau makes the film
    h vary along z (thin end to thick end). Symmetric 4 MPa axial-end BCs. The
    circumferential x-y film is in journal_film_coordinates.png.
    """

    fig, ax = plt.subplots(figsize=(8.6, 6.2), dpi=140)
    fig.suptitle("저널 필름 축 단면 — 틸트 (2D Reynolds 확장, §4.10)", fontsize=12, color="none")

    tau = 0.28
    zc = [-1.0, 1.0]
    yc = [-tau * z for z in zc]
    hh = 0.55
    ax.plot([-1, 1], [1, 1], "--", color="#8a8f98", lw=1.6)
    ax.plot([-1, 1], [-1, -1], "--", color="#8a8f98", lw=1.6)
    ax.fill_between(zc, [y - hh for y in yc], [y + hh for y in yc], color="#dbe7f5", zorder=0)
    ax.plot(zc, [y + hh for y in yc], color="#2166ac", lw=2.2)
    ax.plot(zc, [y - hh for y in yc], color="#2166ac", lw=2.2)
    ax.plot(zc, yc, ":", color="#2166ac", lw=1.4)

    for zz in (-1, 1):
        ax.annotate(
            "",
            xy=(zz, 1.0),
            xytext=(zz, (-tau * zz) + hh),
            arrowprops={"arrowstyle": "<|-|>", "color": "#7a4fb0", "lw": 1.5},
        )
    ax.annotate(r"$h$ (얇음)", (-1, 1), (-1.05, 1.09), color="#7a4fb0", fontsize=10)
    ax.annotate(r"$h$ (두꺼움)", (1, 1), (0.52, 1.09), color="#7a4fb0", fontsize=10)

    ax.plot([0, 0.6], [0, 0], color="#999", lw=1.0)
    ax.add_patch(
        Arc(
            (0, 0),
            0.9,
            0.9,
            angle=0,
            theta1=-degrees(atan2(tau, 1)),
            theta2=0,
            color="#c0392b",
            lw=1.6,
        )
    )
    ax.annotate(
        r"$\tau$ (틸트/미스얼라인먼트)", (0.45, -0.12), (-0.9, -0.5), color="#c0392b", fontsize=10
    )

    ax.annotate(
        "",
        xy=(1.15, -1.35),
        xytext=(-1.15, -1.35),
        arrowprops={"arrowstyle": "-|>", "color": "black"},
    )
    ax.annotate(r"$\zeta=z\ \in[-L_j/2,\,+L_j/2]$", (1.15, -1.35), (-0.6, -1.63), fontsize=11)
    ax.plot([-1, -1], [-1.28, -1.42], color="black", lw=1)
    ax.plot([1, 1], [-1.28, -1.42], color="black", lw=1)
    ax.annotate(r"$L_j=21$ mm", (0, -1.28), (-0.28, -1.20), fontsize=9, color="#555")

    ax.text(
        -1.34,
        -2.0,
        r"$\vec e(\zeta)=\vec e_0+\zeta\,\vec\tau$,  $p(\pm L_j/2)=4$ MPa,  $p_{cav}=4$ MPa"
        "\n"
        r"$h(\alpha,\zeta)=c_j+[(e_x+\zeta\tau_x)\cos\alpha"
        r"+(e_y+\zeta\tau_y)\sin\alpha]$",
        fontsize=10.5,
        bbox={"boxstyle": "round", "fc": "white", "ec": "#c7ccd4"},
    )
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-2.1, 1.6)

    _save(fig, path)


def render_bush_film_coordinates(geometry: RotaryGeometry, path: Path) -> None:
    """Swing-bush film coordinate frames — curved F_c and flat F_f (PHYSICS.md 4.11).

    Panel A: real placement (rotor groove + vane + the two bush pieces) with the
    two ~10 um films marked. Panel B: curved film F_c (partial journal — bush OD
    vs rotor groove). Panel C: flat film F_f (slider — bush flat vs vane). All
    clearances exaggerated; the two bush pieces share one representative F_c/F_f.
    """

    bush = SwingBush()
    state = prescribed_state(geometry, 0.0)
    gx_m, gy_m = state.cutout_center_m
    gx, gy = gx_m / MM, gy_m / MM
    cut_r = geometry.cutout_radius_m / MM
    half_w = 0.5 * geometry.vane_width_m / MM

    fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(16.5, 6.2), dpi=140)
    fig.suptitle("스윙 부시 필름 좌표 — 곡면 $F_c$ · 평면 $F_f$ (§4.11)", fontsize=13, color="none")

    # ---- Panel A: real placement (zoom on the groove) ----
    ax_a.set_aspect("equal")
    ax_a.axis("off")
    ax_a.set_xlim(gx - 13, gx + 13)
    ax_a.set_ylim(gy - 11, gy + 14)
    ax_a.add_patch(
        plt.Circle((gx, gy), cut_r, facecolor="#eef1f5", edgecolor="#888", lw=1.4, zorder=1)
    )
    ax_a.add_patch(
        Rectangle(
            (-half_w, gy - 10.5),
            2 * half_w,
            23.0,
            facecolor=VANE_COLOR,
            edgecolor="black",
            lw=1.0,
            zorder=2,
        )
    )
    for side in (1.0, -1.0):
        xs, ys = _bush_outline_mm(gx_m + side * bush.piece_shift_m, gy_m, side)
        ax_a.fill(xs, ys, facecolor="#9fb6d4", edgecolor="black", lw=1.0, zorder=3)
    ax_a.annotate(
        "평면 필름 $F_f$\n(부시 평면 ↔ 베인) → C",
        xy=(half_w, gy - 2.5),
        xytext=(gx + 4.0, gy - 8.5),
        fontsize=9,
        color="#c0392b",
        arrowprops={"arrowstyle": "->", "color": "#c0392b"},
    )
    ax_a.annotate(
        "곡면 필름 $F_c$\n(부시 OD ↔ 그루브) → B",
        xy=(gx - cut_r + 0.3, gy + 2.0),
        xytext=(gx - 12.5, gy + 9.5),
        fontsize=9,
        color="#2e7d32",
        arrowprops={"arrowstyle": "->", "color": "#2e7d32"},
    )
    ax_a.text(
        gx,
        gy + 12.7,
        "A. 실제 배치 (유막 10 " + r"$\mu$m — 육안 미세)",
        ha="center",
        fontsize=11,
        fontweight="bold",
    )

    # ---- Panel B: curved film F_c (partial journal) ----
    ax_b.set_aspect("equal")
    ax_b.axis("off")
    ax_b.set_xlim(-1.5, 1.5)
    ax_b.set_ylim(-1.7, 1.5)
    r_g, r_p = 1.0, 0.80
    ec = 0.5 * (r_g - r_p)
    eang = radians(78)
    ex, ey = ec * cos(eang), ec * sin(eang)  # e_c = O_p - O_g
    ring = [(cos(radians(t)), sin(radians(t))) for t in range(0, 361, 3)]
    ax_b.plot([r_g * x for x, _ in ring], [r_g * y for _, y in ring], "--", color="#8a8f98", lw=1.4)
    b0, dh = 90, 44
    arc = [
        (ex + r_p * cos(radians(t)), ey + r_p * sin(radians(t)))
        for t in range(b0 - dh, b0 + dh + 1)
    ]
    ax_b.fill(
        [ex] + [p[0] for p in arc],
        [ey] + [p[1] for p in arc],
        facecolor="#dbe7f5",
        edgecolor="none",
        zorder=0,
    )
    ax_b.plot([p[0] for p in arc], [p[1] for p in arc], color="#2166ac", lw=3.0)
    ax_b.annotate("부시 OD (곡면=저널)", (ex, 0.55), (0.28, 0.30), fontsize=9, color="#2166ac")
    ax_b.plot(0, 0, "o", color="#c0392b", ms=6)
    ax_b.annotate(r"$O_g$ (그루브 중심)", (0, 0), (0.08, -0.16), fontsize=9.5, color="#c0392b")
    ax_b.plot(ex, ey, "o", color="#2166ac", ms=5)
    ax_b.annotate(
        "",
        xy=(ex, ey),
        xytext=(0, 0),
        arrowprops={"arrowstyle": "-|>", "color": "#c0392b", "lw": 1.8},
    )
    ax_b.annotate(
        r"$\vec e_c=O_p-O_g$", (ex / 2, ey / 2), (-1.42, 0.02), fontsize=10, color="#c0392b"
    )
    bs = radians(70)
    ax_b.plot([0, r_g * cos(bs)], [0, r_g * sin(bs)], color="#2e7d32", lw=1.1)
    ax_b.plot([0, 0], [0, r_g], color="#2e7d32", lw=0.8, ls=":")
    ax_b.add_patch(Arc((0, 0), 1.1, 1.1, angle=0, theta1=70, theta2=90, color="#2e7d32", lw=1.5))
    ax_b.annotate(r"$\beta$", (0.17, 0.62), (0.17, 0.62), color="#2e7d32", fontsize=12)
    gx0, gy0 = r_g * cos(radians(98)), r_g * sin(radians(98))
    bx0, by0 = ex + r_p * cos(radians(98)), ey + r_p * sin(radians(98))
    ax_b.annotate(
        "",
        xy=(gx0, gy0),
        xytext=(bx0, by0),
        arrowprops={"arrowstyle": "<|-|>", "color": "#7a4fb0", "lw": 1.6},
    )
    ax_b.annotate(
        r"$h_c(\beta)$", (gx0, gy0), (gx0 - 0.30, gy0 + 0.06), color="#7a4fb0", fontsize=11
    )
    ax_b.annotate(
        "",
        xy=(r_g * cos(radians(122)), r_g * sin(radians(122))),
        xytext=(r_g * cos(radians(107)), r_g * sin(radians(107))),
        arrowprops={"arrowstyle": "-|>", "color": "black", "lw": 1.6},
    )
    ax_b.annotate(
        "$U_c$ (로터 스윙)",
        (r_g * cos(radians(130)), r_g * sin(radians(130))),
        (-1.45, 1.02),
        fontsize=9,
    )
    ax_b.text(0.0, -0.75, "로터 그루브 (오목=보어)", ha="center", fontsize=9, color="#666")
    ax_b.text(-1.45, 1.34, "B. 곡면 $F_c$ — 부분 저널", fontsize=11, fontweight="bold")
    ax_b.text(
        -1.48,
        -1.64,
        r"$h_c(\beta)=c_c-\vec e_c\!\cdot\!\hat r$,  $c_c\approx10\,\mu$m,  축 $\zeta$"
        "\n"
        r"BC: 챔버압 ↔ $p_{rec}=4$ MPa;  전단=로터 스윙",
        fontsize=9,
        bbox={"boxstyle": "round", "fc": "white", "ec": "#c7ccd4"},
    )

    # ---- Panel C: flat film F_f (slider) ----
    ax_c.set_aspect("equal")
    ax_c.axis("off")
    ax_c.set_xlim(-1.5, 2.3)
    ax_c.set_ylim(-1.7, 1.6)
    cf = 0.16
    ax_c.add_patch(
        Rectangle((-1.2, -1.0), 1.2, 2.2, facecolor=VANE_COLOR, edgecolor="black", lw=1.0)
    )
    ax_c.text(-0.6, 1.30, "베인 (고정)", ha="center", fontsize=10)
    ax_c.add_patch(
        Rectangle((cf, -0.85), 1.1, 1.9, facecolor="#dbe7f5", edgecolor="#2166ac", lw=2.0)
    )
    ax_c.text(cf + 0.55, 0.1, "부시\n평면", ha="center", fontsize=10, color="#2166ac")
    ax_c.add_patch(
        Rectangle((0.0, -0.85), cf, 1.9, facecolor="#fbeaea", edgecolor="none", zorder=0)
    )
    ax_c.annotate(
        "",
        xy=(cf, -0.55),
        xytext=(0.0, -0.55),
        arrowprops={"arrowstyle": "<|-|>", "color": "#7a4fb0", "lw": 1.6},
    )
    ax_c.annotate(r"$h_f(s)$", (cf, -0.55), (cf + 0.05, -0.64), color="#7a4fb0", fontsize=11)
    ax_c.annotate(
        "",
        xy=(0.0, 1.15),
        xytext=(0.0, -0.95),
        arrowprops={"arrowstyle": "-|>", "color": "#2e7d32", "lw": 1.4},
    )
    ax_c.annotate(
        r"$s$ (베인 따라, $0..L_f$)", (0.0, 1.15), (0.06, 1.22), color="#2e7d32", fontsize=9.5
    )
    ax_c.annotate(
        "",
        xy=(cf + 0.55, 1.25),
        xytext=(cf + 0.55, 0.55),
        arrowprops={"arrowstyle": "-|>", "color": "black", "lw": 2.0},
    )
    ax_c.annotate("$U_f$ (부시 병진)", (cf + 0.55, 1.25), (cf + 0.72, 0.8), fontsize=9)
    ax_c.plot(1.7, -1.3, marker="o", mfc="white", mec="black", ms=13)
    ax_c.plot(1.7, -1.3, marker=".", color="black", ms=4)
    ax_c.annotate(r"$\zeta=z\odot$ (축)", (1.7, -1.3), (0.35, -1.35), fontsize=9.5)
    ax_c.text(-1.45, 1.48, "C. 평면 $F_f$ — 슬라이더", fontsize=11, fontweight="bold")
    ax_c.text(
        -1.48,
        -1.64,
        r"$h_f(s)=c_f-\Delta x+s\,\gamma$ (틸트 시 쐐기),  $c_f\approx10\,\mu$m"
        "\n"
        r"BC: 챔버압 ↔ $p_{rec}=4$ MPa;  전단=부시 병진",
        fontsize=9,
        bbox={"boxstyle": "round", "fc": "white", "ec": "#c7ccd4"},
    )

    _save(fig, path)


def render_bush_attitude_reference(geometry: RotaryGeometry, path: Path) -> None:
    """Vane-referenced attitude angles for the moving swing bush (PHYSICS.md 4.11).

    Panel A: real placement at max swing (theta=90) — the rotor attitude phi_r is
    the O_r->O_g axis measured from the vane (ground), not the bush. Panel B: the
    angle relations (exaggerated) — vane datum, rotor phi_r, bush phi_b, and the
    curved-film relative swing phi_r - phi_b. Reduces to Section 3.6/4.7 at phi_b=0.
    """

    theta = radians(90.0)
    state = prescribed_state(geometry, theta)
    orient = state.rotor_orientation_rad
    rx, ry = state.rotor_center_m[0] / MM, state.rotor_center_m[1] / MM
    gx, gy = state.cutout_center_m[0] / MM, state.cutout_center_m[1] / MM
    bore = geometry.cylinder_radius_m / MM
    half_w = 0.5 * geometry.vane_width_m / MM
    vtip = state.vane_tip_m[1] / MM
    lim = bore + 5

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14.0, 7.0), dpi=140)
    fig.suptitle("베인 기준 자세각 — 이동하는 스윙 부시 (§4.11)", fontsize=13, color="none")

    # ---- Panel A: real placement at max swing ----
    ax_a.set_aspect("equal")
    ax_a.axis("off")
    ax_a.set_xlim(-lim, lim)
    ax_a.set_ylim(-lim, lim)
    ax_a.add_patch(
        plt.Circle((0, 0), bore, facecolor="#f7f8fa", edgecolor="black", lw=1.6, zorder=1)
    )
    contour = rotor_contour(geometry, theta)
    ax_a.fill(
        [p[0] / MM for p in contour.material],
        [p[1] / MM for p in contour.material],
        facecolor=ROTOR_COLOR,
        edgecolor="black",
        lw=1.3,
        zorder=2,
    )
    bush = SwingBush()
    for side in (1.0, -1.0):
        xs, ys = _bush_outline_mm(
            state.cutout_center_m[0] + side * bush.piece_shift_m, state.cutout_center_m[1], side
        )
        ax_a.fill(xs, ys, facecolor=BUSH_COLOR, edgecolor="black", lw=0.8, zorder=3)
    ax_a.add_patch(
        Rectangle(
            (-half_w, vtip),
            2 * half_w,
            bore - vtip,
            facecolor=VANE_COLOR,
            edgecolor="black",
            lw=1.0,
            zorder=4,
        )
    )
    ax_a.plot(0, 0, "+", color="black", ms=11, mew=1.8, zorder=6)
    ax_a.plot(rx, ry, "o", color="#c0392b", ms=5, zorder=7)
    ax_a.annotate(r"$O_r$", (rx, ry), (rx + 1.3, ry - 1.6), fontsize=10)
    ax_a.plot(gx, gy, "o", color="#2166ac", ms=5, zorder=7)
    ax_a.annotate(r"$O_g$", (gx, gy), (gx + 1.0, gy + 0.4), fontsize=10, color="#2166ac")
    ax_a.plot([rx, rx], [ry, ry + 13], "--", color="#888", lw=1.3, zorder=5)
    ax_a.annotate(
        r"베인 방향 ($+\hat y$)", (rx, ry + 13), (rx - 13, ry + 12.2), fontsize=9, color="#888"
    )
    axis_len = 13
    ax_a.annotate(
        "",
        xy=(rx + axis_len * cos(orient), ry + axis_len * sin(orient)),
        xytext=(rx, ry),
        arrowprops={"arrowstyle": "-|>", "color": "#c0392b", "lw": 1.9},
        zorder=6,
    )
    ax_a.annotate(
        r"로터 자세축 $O_r\!\to\!O_g$", (rx, ry), (rx + 2.0, ry + 11.0), fontsize=9, color="#c0392b"
    )
    ax_a.add_patch(
        Arc((rx, ry), 15, 15, angle=0, theta1=90, theta2=degrees(orient), color="#7a4fb0", lw=2.0)
    )
    ax_a.annotate(
        rf"$\phi_r={degrees(orient) - 90:+.1f}°$",
        (rx, ry),
        (rx - 6.5, ry + 6.0),
        fontsize=12,
        color="#7a4fb0",
    )
    ax_a.text(
        rx - 1,
        ry - 9,
        r"부시는 베인과 정렬 ($\phi_b\approx0$, 명목)",
        ha="center",
        fontsize=8.5,
        color="#555",
    )
    ax_a.text(
        -lim + 1,
        lim - 2,
        "A. 실제 배치 (θ=90°, 최대 스윙) — 로터 자세 $\\phi_r$",
        fontsize=11.5,
        fontweight="bold",
    )

    # ---- Panel B: angle relations (exaggerated) ----
    ax_b.set_aspect("equal")
    ax_b.axis("off")
    ax_b.set_xlim(-1.5, 1.5)
    ax_b.set_ylim(-1.6, 1.35)
    phr, phb = radians(30), radians(11)

    def _ray(ang: float, r: float) -> tuple[float, float]:
        return (-r * sin(ang), r * cos(ang))  # tilt left (CCW) from +y

    ax_b.annotate(
        "",
        xy=(0, 1.15),
        xytext=(0, 0),
        arrowprops={"arrowstyle": "-|>", "color": "#888", "lw": 1.9},
    )
    ax_b.annotate("베인 (기준)", (0, 1.15), (0.04, 1.2), fontsize=10, color="#555")
    ax_b.annotate(
        "",
        xy=_ray(phr, 1.15),
        xytext=(0, 0),
        arrowprops={"arrowstyle": "-|>", "color": "#c0392b", "lw": 2.1},
    )
    ax_b.annotate(
        r"로터 자세 $\phi_r$",
        _ray(phr, 1.15),
        (_ray(phr, 1.15)[0] - 0.55, _ray(phr, 1.15)[1] + 0.02),
        fontsize=10,
        color="#c0392b",
    )
    ax_b.annotate(
        "",
        xy=_ray(phb, 1.15),
        xytext=(0, 0),
        arrowprops={"arrowstyle": "-|>", "color": "#2166ac", "lw": 2.1},
    )
    ax_b.annotate(
        r"부시 자세 $\phi_b$",
        _ray(phb, 1.15),
        (_ray(phb, 1.15)[0] + 0.06, _ray(phb, 1.15)[1] - 0.05),
        fontsize=10,
        color="#2166ac",
    )
    ax_b.add_patch(
        Arc((0, 0), 0.7, 0.7, angle=0, theta1=90, theta2=90 + degrees(phb), color="#2166ac", lw=1.7)
    )
    ax_b.annotate(
        r"$\phi_b$",
        _ray(phb / 2, 0.42),
        (_ray(phb / 2, 0.42)[0] - 0.02, _ray(phb / 2, 0.42)[1]),
        fontsize=11,
        color="#2166ac",
    )
    ax_b.add_patch(
        Arc((0, 0), 1.4, 1.4, angle=0, theta1=90, theta2=90 + degrees(phr), color="#c0392b", lw=1.7)
    )
    ax_b.annotate(
        r"$\phi_r$",
        _ray(phr * 0.55, 0.78),
        (_ray(phr * 0.55, 0.78)[0] - 0.14, _ray(phr * 0.55, 0.78)[1]),
        fontsize=11,
        color="#c0392b",
    )
    ax_b.add_patch(
        Arc(
            (0, 0),
            1.9,
            1.9,
            angle=0,
            theta1=90 + degrees(phb),
            theta2=90 + degrees(phr),
            color="#2e7d32",
            lw=2.2,
        )
    )
    mid = (phb + phr) / 2
    ax_b.annotate(
        r"$\phi_r-\phi_b$" + "\n(곡면 상대 스윙)",
        _ray(mid, 1.15),
        (_ray(mid, 1.15)[0] - 0.15, _ray(mid, 1.15)[1] + 0.16),
        fontsize=9.5,
        color="#2e7d32",
        ha="center",
    )
    ax_b.text(-1.48, 1.28, "B. 각 관계 (과장)", fontsize=11.5, fontweight="bold")
    ax_b.text(
        -1.5,
        -1.55,
        r"$\phi_r=\phi-90°$ (베인 기준 로터 스윙 $\pm10.37°$);  "
        r"$\phi_b$ 부시 자세 (0 → DOF)"
        "\n"
        r"곡면 $F_c$: $U_c=r_b(\dot\phi_r-\dot\phi_b)$;  "
        r"평면 $F_f$: 쐐기 $\gamma=\phi_b$;  환원 $\phi_b\!=\!0\Rightarrow$ §3.6/4.7",
        fontsize=9,
        bbox={"boxstyle": "round", "fc": "white", "ec": "#c7ccd4"},
    )

    _save(fig, path)


def render_reynolds_1d_validation(geometry: RotaryGeometry, path: Path) -> None:
    """1-D numerical Reynolds solver vs the Ocvirk closed form (PHYSICS.md 4.12).

    Panel A: load capacity |F| over the eccentricity ratio — the Ocvirk closed form
    (line) and the 1-D finite-difference Reynolds solve (points) coincide. Panel B:
    the relative error, ~1e-4 (static; the residual is only the circumferential
    quadrature since the axial solution is an exact parabola).
    """

    omega = 95.0
    eps = [0.05 * i for i in range(1, 19)]  # 0.05 .. 0.90
    ocv = [short_bearing_force(e, omega).magnitude_n for e in eps]
    num = [solve_short_bearing_1d(e, omega).magnitude_n for e in eps]
    rel = [abs(n - o) / o for n, o in zip(num, ocv, strict=True)]

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13.2, 5.4), dpi=140)
    fig.suptitle("1D 수치 Reynolds 솔버 검증 — Ocvirk 해석해 (§4.12)", fontsize=12, color="none")

    ax_a.semilogy(eps, ocv, color="#2166ac", lw=2.2, label="Ocvirk 닫힌 형태 (§4.9)")
    ax_a.semilogy(
        eps, num, "o", color="#c0392b", ms=6, mfc="none", mew=1.6, label="1D 수치 Reynolds (§4.12)"
    )
    ax_a.set_xlabel("편심비 ε")
    ax_a.set_ylabel("하중 용량 |F| (N)")
    ax_a.grid(True, which="both", color="#e2e5ea", lw=0.7)
    ax_a.text(
        0.02,
        1.02,
        "A. 하중 용량 (정적, 순수 회전)",
        transform=ax_a.transAxes,
        fontsize=10,
        fontweight="bold",
    )
    _legend_above(ax_a, ncol=1)

    ax_b.plot(eps, [r * 1.0e4 for r in rel], "o-", color="#2e7d32", lw=1.6, ms=5)
    ax_b.set_xlabel("편심비 ε")
    ax_b.set_ylabel("상대 오차 (×1e-4)")
    ax_b.set_ylim(bottom=0.0)
    ax_b.grid(True, color="#e2e5ea", lw=0.7)
    ax_b.text(
        0.02,
        1.02,
        "B. Ocvirk 대비 상대 오차",
        transform=ax_b.transAxes,
        fontsize=10,
        fontweight="bold",
    )
    ax_b.text(
        0.03,
        0.05,
        "정적 ~1e-4 (축 해=포물선, 잔차=원주 적분)\n"
        "스퀴즈 시 ~2%: π-film 중첩 vs 실제 캐비테이션역",
        transform=ax_b.transAxes,
        fontsize=8.5,
        va="bottom",
        bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "#c7ccd4"},
    )
    _save(fig, path)


def render_rotor_orbit(curve: PressureCurve, path: Path) -> None:
    """Rotor lateral orbit — dynamic EOM vs quasi-static (PHYSICS.md 4.13).

    Panel A: eccentricity ratio over the crank angle — the dynamic time-integrated
    orbit (§4.13) vs the quasi-static estimate (§4.9); the squeeze-film lag
    attenuates and phase-lags the swing (peak ~0.50 vs ~0.71). Panel B: the whirl
    orbit (e_x, e_y) inside the journal clearance circle.
    """

    orbit = integrate_rotor_orbit(curve.geometry, trace=curve.trace)
    ang = [degrees(a) for a in orbit.crank_angle_rad]
    order = sorted(range(len(ang)), key=lambda i: ang[i])
    theta = [ang[i] for i in order]
    eps_dyn = [orbit.eccentricity_ratio[i] for i in order]
    eps_qs = [orbit.quasi_static_eccentricity_ratio[i] for i in order]
    clr_um = orbit.clearance_m * 1e6
    ex_um = [e * 1e6 for e in orbit.eccentricity_x_m]
    ey_um = [e * 1e6 for e in orbit.eccentricity_y_m]

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13.4, 5.6), dpi=140)
    fig.suptitle("로터 횡방향 궤도 — 동역학 EOM vs 준정적 (§4.13)", fontsize=12, color="none")

    ax_a.plot(
        theta, eps_qs, "--", color="#9aa0a8", lw=2.0, label=r"준정적 (§4.9, $\dot\varepsilon$=0)"
    )
    ax_a.plot(theta, eps_dyn, "-", color="#c0392b", lw=2.3, label="동역학 EOM (§4.13)")
    ax_a.set_xlabel("크랭크 각 θ (deg)")
    ax_a.set_ylabel("편심비 ε")
    ax_a.set_xlim(0, 360)
    ax_a.set_xticks(range(0, 361, 90))
    ax_a.set_ylim(0, 1)
    ax_a.grid(True, color="#e2e5ea", lw=0.7)
    ax_a.text(
        0.02,
        1.02,
        "A. 편심비 ε(θ) — 스퀴즈 감쇠 지연",
        transform=ax_a.transAxes,
        fontsize=10,
        fontweight="bold",
    )
    ax_a.text(
        0.03,
        0.05,
        f"동역학 피크 {orbit.peak_eccentricity_ratio:.2f} "
        f"< 준정적 {orbit.quasi_static_peak_eccentricity_ratio:.2f}\n"
        f"최소 유막 {orbit.minimum_film_thickness_m * 1e6:.1f} µm; 원심 "
        f"{orbit.centrifugal_load_n:.0f} N (~2%)\n"
        f"저널 마찰: 동역학 {orbit.dynamic_journal_friction_w:.1f}W "
        f"(Petroff {orbit.petroff_journal_friction_w:.1f} / 준정적 "
        f"{orbit.quasi_static_journal_friction_w:.1f})",
        transform=ax_a.transAxes,
        fontsize=8.5,
        va="bottom",
        bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "#c7ccd4"},
    )
    _legend_above(ax_a, ncol=2)

    ax_b.set_aspect("equal")
    ring = [(cos(radians(t)), sin(radians(t))) for t in range(0, 361, 2)]
    ax_b.plot(
        [clr_um * x for x, _ in ring],
        [clr_um * y for _, y in ring],
        "--",
        color="#8a8f98",
        lw=1.4,
        label=f"클리어런스 원 (c_j={clr_um:.0f} µm)",
    )
    ax_b.plot(
        ex_um + [ex_um[0]],
        ey_um + [ey_um[0]],
        "-",
        color="#c0392b",
        lw=2.0,
        label="로터 whirl 궤도",
    )
    ax_b.plot(0, 0, "+", color="black", ms=10, mew=1.6)
    ax_b.set_xlabel(r"$e_x$ (µm)")
    ax_b.set_ylabel(r"$e_y$ (µm)")
    ax_b.grid(True, color="#e2e5ea", lw=0.7)
    ax_b.text(
        0.02,
        1.02,
        "B. whirl 궤도 (클리어런스 내)",
        transform=ax_b.transAxes,
        fontsize=10,
        fontweight="bold",
    )
    ax_b.legend(loc="lower center", bbox_to_anchor=(0.5, 1.09), ncol=1, fontsize=8, framealpha=0.9)
    _save(fig, path)


def render_bush_2piece_coordinates(geometry: RotaryGeometry, path: Path) -> None:
    """Two independent swing-bush pieces — DOF and film frames (PHYSICS.md 4.11).

    Panel A: real placement — the IN and OUT pieces flank the vane in the rotor
    groove, each an independent body with two free DOF (vane-normal u_k, attitude
    phi_k) and two films (curved F_c vs groove wall, flat F_f vs vane). Panel B: a
    single piece's film frames (curved arc beta, flat slider s). Gaps exaggerated.
    """

    bush = SwingBush()
    state = prescribed_state(geometry, 0.0)
    gx_m, gy_m = state.cutout_center_m
    gx, gy = gx_m / MM, gy_m / MM
    cut_r = geometry.cutout_radius_m / MM
    half_w = 0.5 * geometry.vane_width_m / MM

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14.4, 7.2), dpi=140)
    fig.suptitle("스윙 부시 — 두 독립 피스 자유도·필름 좌표 (§4.11)", fontsize=13, color="none")

    # ---- Panel A: real placement, per-piece DOF ----
    ax_a.set_aspect("equal")
    ax_a.axis("off")
    ax_a.set_xlim(gx - 15, gx + 15)
    ax_a.set_ylim(gy - 13, gy + 13)
    ax_a.add_patch(
        plt.Circle((gx, gy), cut_r, facecolor="#eef1f5", edgecolor="#888", lw=1.4, zorder=1)
    )
    ax_a.add_patch(
        Rectangle(
            (-half_w, gy - 12),
            2 * half_w,
            26,
            facecolor=VANE_COLOR,
            edgecolor="black",
            lw=1.0,
            zorder=2,
        )
    )
    ax_a.text(0, gy + 11.6, "베인 (고정)", ha="center", fontsize=9, color="white")
    ax_a.text(
        gx,
        gy - 11.2,
        "로터 그루브 (반경 8mm) — 로터와 함께 이동",
        ha="center",
        fontsize=8.5,
        color="#666",
    )
    ax_a.plot(gx, gy, "o", color="white", mec="#c0392b", mew=1.6, ms=8, zorder=6)
    ax_a.annotate(r"$O_g$", (gx, gy), (gx + 0.6, gy + 0.6), fontsize=10, color="#c0392b", zorder=7)
    for side, tag, col, fill in (
        (1.0, "IN", "#2166ac", "#cdd9ea"),
        (-1.0, "OUT", "#2e7d32", "#cfe3d3"),
    ):
        xs, ys = _bush_outline_mm(gx_m + side * bush.piece_shift_m, gy_m, side)
        ax_a.fill(xs, ys, facecolor=fill, edgecolor="black", lw=1.1, zorder=3)
        pcx = (
            gx
            + side
            * (bush.flat_offset_m + 0.45 * (bush.piece_outer_radius_m - bush.flat_offset_m))
            / MM
        )
        pcy = gy
        ax_a.text(
            pcx, pcy + 3.2, f"{tag} 피스", ha="center", fontsize=9, color=col, fontweight="bold"
        )
        ax_a.annotate(
            "",
            xy=(pcx + side * 2.6, pcy - 0.3),
            xytext=(pcx - side * 0.6, pcy - 0.3),
            arrowprops={"arrowstyle": "<|-|>", "color": col, "lw": 1.8},
            zorder=6,
        )
        ax_a.annotate(
            rf"$u_{{{tag}}}$",
            (pcx + side, pcy - 0.3),
            (pcx + side * 0.9, pcy - 2.2),
            fontsize=9.5,
            color=col,
            ha="center",
        )
        ax_a.add_patch(
            Arc((pcx, pcy), 3.6, 3.6, angle=0, theta1=30, theta2=95, color=col, lw=1.7, zorder=6)
        )
        ax_a.annotate(
            rf"$\phi_{{{tag}}}$",
            (pcx, pcy + 1.9),
            (pcx - side * 1.4, pcy + 1.6),
            fontsize=10,
            color=col,
            ha="center",
        )
        ax_a.annotate(
            rf"곡면 $F_c$ ({tag})",
            (gx + side * (cut_r - 0.3), gy + 2.0),
            (gx + side * 10.5, gy + 7.5),
            fontsize=8.5,
            color=col,
            ha="center",
            arrowprops={"arrowstyle": "->", "color": col, "lw": 1.1},
        )
        ax_a.annotate(
            rf"평면 $F_f$ ({tag})",
            (side * half_w, gy - 4.0),
            (gx + side * 9.5, gy - 8.0),
            fontsize=8.5,
            color=col,
            ha="center",
            arrowprops={"arrowstyle": "->", "color": col, "lw": 1.1},
        )
    ax_a.text(
        gx - 14.5,
        gy + 12.0,
        r"A. 두 피스 독립 — 각 ($u_k$, $\phi_k$) 2 DOF",
        fontsize=11,
        fontweight="bold",
    )

    # ---- Panel B: single piece (IN) film frames ----
    ax_b.set_aspect("equal")
    ax_b.axis("off")
    ax_b.set_xlim(-2.6, 3.0)
    ax_b.set_ylim(-2.6, 2.4)
    r_p, cf, cc = 1.7, 0.32, 0.34
    arc = [radians(a) for a in range(-90, 91, 3)]
    ax_b.fill(
        [0.0] + [r_p * cos(t) for t in arc] + [0.0],
        [r_p] + [r_p * sin(t) for t in arc] + [-r_p],
        facecolor="#cdd9ea",
        edgecolor="#2166ac",
        lw=2.2,
        zorder=1,
    )
    ax_b.text(0.65, -0.95, "IN 피스", ha="center", fontsize=9.5, color="#2166ac", fontweight="bold")
    ax_b.add_patch(
        Rectangle((-cf - 0.9, -1.9), 0.9, 3.8, facecolor=VANE_COLOR, edgecolor="black", lw=1.0)
    )
    ax_b.text(-cf - 0.45, 2.05, "베인", ha="center", fontsize=9)
    ax_b.add_patch(Rectangle((-cf, -r_p), cf, 2 * r_p, facecolor="#fbeaea", zorder=0))
    ax_b.annotate(
        "",
        xy=(0, -1.15),
        xytext=(-cf, -1.15),
        arrowprops={"arrowstyle": "<|-|>", "color": "#7a4fb0", "lw": 1.5},
    )
    ax_b.annotate(
        r"$h_f$",
        (-cf / 2, -1.15),
        (-cf / 2 - 0.05, -1.55),
        fontsize=10,
        color="#7a4fb0",
        ha="center",
    )
    ax_b.annotate(
        "",
        xy=(-cf - 0.02, 1.7),
        xytext=(-cf - 0.02, -1.7),
        arrowprops={"arrowstyle": "-|>", "color": "#7a4fb0", "lw": 1.2},
    )
    ax_b.annotate(
        r"평면 $F_f$: $s$", (-cf - 0.02, 1.7), (-cf - 0.95, 1.72), fontsize=9, color="#7a4fb0"
    )
    wall = [radians(a) for a in range(-84, 85, 3)]
    ax_b.plot(
        [(r_p + cc) * cos(t) for t in wall],
        [(r_p + cc) * sin(t) for t in wall],
        "--",
        color="#8a8f98",
        lw=1.6,
    )
    ax_b.text(r_p + cc + 0.12, 1.15, "그루브 벽", ha="left", fontsize=9, color="#666")
    ax_b.annotate(
        "",
        xy=((r_p + cc) * cos(radians(28)), (r_p + cc) * sin(radians(28))),
        xytext=(r_p * cos(radians(28)), r_p * sin(radians(28))),
        arrowprops={"arrowstyle": "<|-|>", "color": "#c0392b", "lw": 1.5},
    )
    ax_b.annotate(r"$h_c$", (r_p + 0.15, 1.05), (r_p + 0.15, 1.05), fontsize=10, color="#c0392b")
    ax_b.add_patch(Arc((0, 0), 1.15, 1.15, angle=0, theta1=-24, theta2=24, color="#c0392b", lw=1.5))
    ax_b.annotate(r"곡면 $F_c$: $\beta$", (0.62, 0.0), (0.9, 0.55), fontsize=9, color="#c0392b")
    ax_b.annotate(
        "",
        xy=(1.1, -2.15),
        xytext=(-0.3, -2.15),
        arrowprops={"arrowstyle": "<|-|>", "color": "black", "lw": 1.7},
    )
    ax_b.annotate(r"$u_{IN}$ (법선)", (0.4, -2.15), (0.4, -2.5), fontsize=9.5, ha="center")
    ax_b.add_patch(Arc((0, 0), 0.62, 0.62, angle=0, theta1=205, theta2=255, color="black", lw=1.5))
    ax_b.annotate(r"$\phi_{IN}$", (0, 0), (-0.62, -0.45), fontsize=10)
    ax_b.text(
        -2.55, 2.15, "B. 단일 피스 필름 프레임 (곡면 β + 평면 s)", fontsize=11, fontweight="bold"
    )

    _save(fig, path)


def _bush_zoom_base(ax, geometry: RotaryGeometry, theta_deg: float) -> tuple[float, float]:
    """Draw the zoomed rotor-mouth / bush / vane region; return the groove centre (mm)."""

    angle = radians(theta_deg)
    state = prescribed_state(geometry, angle)
    cx = state.cutout_center_m[0] / MM
    cy = state.cutout_center_m[1] / MM
    contour = rotor_contour(geometry, angle)
    ax.fill(
        [p[0] / MM for p in contour.material],
        [p[1] / MM for p in contour.material],
        facecolor=ROTOR_COLOR,
        edgecolor="black",
        linewidth=1.2,
        zorder=1,
    )
    ax.add_patch(
        plt.Circle(
            (cx, cy),
            geometry.cutout_radius_m / MM,
            facecolor="#eef1f5",
            edgecolor="#888",
            linewidth=1.2,
            zorder=1,
        )
    )
    bush = SwingBush()
    for side in (1.0, -1.0):
        xs, ys = _bush_outline_mm(
            state.cutout_center_m[0] + side * bush.piece_shift_m, state.cutout_center_m[1], side
        )
        ax.fill(xs, ys, facecolor=BUSH_COLOR, edgecolor="black", linewidth=0.9, zorder=3)
    vane_x, vane_y = _vane_outline_mm(geometry, state.vane_tip_m[1])
    ax.fill(vane_x, vane_y, facecolor=VANE_COLOR, edgecolor="black", linewidth=1.0, zorder=2)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(cx - 13.0, cx + 13.0)
    ax.set_ylim(cy - 13.0, cy + 12.0)
    return cx, cy


def render_rotor_swing_vs_bush(geometry: RotaryGeometry, path: Path) -> None:
    """Zoom defining the rotor's angular velocity relative to the swing bush (§3.6, 곡면 필름).

    The bush keeps a fixed attitude and only the rotor swings, so the curved film
    (bush curved face <-> rotor groove) is sheared by the rotor-vs-bush relative
    rotation ``omega_rel = omega * dphi/dtheta``.
    """

    theta_deg = 45.0
    fig, ax = plt.subplots(figsize=(7.4, 7.4), dpi=140)
    cx, cy = _bush_zoom_base(ax, geometry, theta_deg)
    fig.suptitle("스윙 부시에 대한 로터 상대 각속도 — 곡면 유막 (§3.6)", color="none")
    groove_r = geometry.cutout_radius_m / MM

    # Rotor swings relative to the fixed-attitude bush: curved double arrow.
    rr = groove_r + 3.0
    ax.annotate(
        "",
        xy=(cx + rr * cos(radians(-35)), cy + rr * sin(radians(-35))),
        xytext=(cx + rr * cos(radians(35)), cy + rr * sin(radians(35))),
        arrowprops={
            "arrowstyle": "<|-|>",
            "color": "#c0392b",
            "lw": 2.2,
            "connectionstyle": "arc3,rad=0.38",
        },
        zorder=6,
    )
    ax.text(cx + rr + 0.8, cy, "로터 스윙", color="#c0392b", fontsize=11, va="center", zorder=7)
    # Curved-film shear (tangential slide) at the bush-outer <-> groove-wall interface.
    fx = cx - groove_r + 0.4
    ax.annotate(
        "",
        xy=(fx, cy - 4.5),
        xytext=(fx, cy + 1.5),
        arrowprops={"arrowstyle": "-|>", "color": "#7a4fb0", "lw": 2.0},
        zorder=6,
    )
    ax.text(fx - 0.6, cy - 6.0, "곡면 유막", color="#7a4fb0", fontsize=10, ha="center", zorder=7)
    ax.text(
        cx - groove_r - 1.5, cy + 6.5, "부시", color="#274060", fontsize=10, ha="right", zorder=7
    )
    _save(fig, path)


def render_bush_translation_vs_vane(geometry: RotaryGeometry, path: Path) -> None:
    """Zoom defining the bush translation speed relative to the vane (§3.6, 평면 필름).

    The vane is fixed; the bush centre follows the groove centre and slides along
    the vane flank, shearing the flat film at ``U_flat = omega * d(y_groove)/dtheta``.
    """

    theta_deg = 45.0
    fig, ax = plt.subplots(figsize=(7.4, 7.4), dpi=140)
    cx, cy = _bush_zoom_base(ax, geometry, theta_deg)
    fig.suptitle("베인에 대한 병진 운동 속도 — 평면 유막 (§3.6)", color="none")
    half_w = 0.5 * geometry.vane_width_m / MM

    # Bush translates along the vane (vertical): double arrow beside the flat face.
    ax.annotate(
        "",
        xy=(half_w + 3.2, cy + 6.0),
        xytext=(half_w + 3.2, cy - 6.0),
        arrowprops={"arrowstyle": "<|-|>", "color": "#c0392b", "lw": 2.4},
        zorder=6,
    )
    ax.text(half_w + 4.0, cy, "부시 병진", color="#c0392b", fontsize=11, va="center", zorder=7)
    # Flat film location (bush flat <-> vane flank).
    ax.annotate(
        "평면 유막",
        xy=(half_w, cy - 1.5),
        xytext=(-half_w - 8.5, cy - 8.0),
        color="#7a4fb0",
        fontsize=10,
        zorder=7,
        arrowprops={"arrowstyle": "->", "color": "#7a4fb0"},
    )
    ax.text(0.0, cy + 10.5, "베인", color="#274060", fontsize=10, ha="center", zorder=7)
    _save(fig, path)


def render_roller_vs_crank(geometry: RotaryGeometry, path: Path) -> None:
    """Schematic defining the roller's angular velocity relative to the crank pin (§4.7).

    The crank pin turns at the shaft speed ``omega`` and the roller (rotor) only
    swings at ``omega*dphi/dtheta``, so the journal film is sheared at their
    difference ``omega_rel = omega(1 - dphi/dtheta)``.
    """

    r_pin = 10.0
    c_draw = 1.1  # exaggerated clearance
    fig, ax = plt.subplots(figsize=(7.4, 7.4), dpi=140)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-16, 16)
    ax.set_ylim(-16, 15)
    fig.suptitle("크랭크 축에 대한 롤러 상대 각속도 — 저널 (§4.7)", color="none")

    # Roller bore ring + oil film + crank pin, concentric.
    ax.add_patch(
        plt.Circle(
            (0, 0),
            r_pin + c_draw + 3.0,
            facecolor=ROTOR_COLOR,
            edgecolor="black",
            linewidth=1.4,
            zorder=1,
        )
    )
    ax.add_patch(
        plt.Circle(
            (0, 0),
            r_pin + c_draw,
            facecolor="#fde8c8",
            edgecolor="#c0392b",
            linewidth=1.2,
            zorder=2,
        )
    )
    ax.add_patch(
        plt.Circle((0, 0), r_pin, facecolor="#6b7a8f", edgecolor="black", linewidth=1.5, zorder=3)
    )
    ax.plot(0, 0, "+", color="white", markersize=9, markeredgewidth=1.6, zorder=4)
    # Crank pin spin omega (fast) inside the pin.
    ax.annotate(
        "",
        xy=(0.55 * r_pin * cos(radians(20)), 0.55 * r_pin * sin(radians(20))),
        xytext=(0.55 * r_pin * cos(radians(160)), 0.55 * r_pin * sin(radians(160))),
        arrowprops={
            "arrowstyle": "-|>",
            "color": "white",
            "lw": 2.2,
            "connectionstyle": "arc3,rad=-0.45",
        },
        zorder=5,
    )
    ax.text(
        0,
        -0.35 * r_pin,
        "크랭크핀",
        color="white",
        fontsize=10,
        ha="center",
        va="center",
        zorder=6,
    )
    # Roller spin (slow, ± swing) on the rotor ring.
    ro = r_pin + c_draw + 3.0
    ax.annotate(
        "",
        xy=(ro * cos(radians(58)), ro * sin(radians(58))),
        xytext=(ro * cos(radians(92)), ro * sin(radians(92))),
        arrowprops={
            "arrowstyle": "-|>",
            "color": "#2e7d32",
            "lw": 2.2,
            "connectionstyle": "arc3,rad=-0.3",
        },
        zorder=6,
    )
    ax.text(
        ro * cos(radians(70)) + 1.0,
        ro * sin(radians(70)) + 1.0,
        "롤러 스핀",
        color="#2e7d32",
        fontsize=10,
        zorder=7,
    )
    ax.annotate(
        "저널 유막",
        xy=(r_pin + 0.5 * c_draw, 0.0),
        xytext=(ro + 0.3, 4.5),
        color="#c0392b",
        fontsize=10,
        ha="left",
        zorder=7,
        arrowprops={"arrowstyle": "->", "color": "#c0392b"},
    )
    _save(fig, path)


# --------------------------------------------------------------------------
# Figure: gap leakage and volumetric efficiency
# --------------------------------------------------------------------------


def render_recompression_pressure(curve: PressureCurve, path: Path) -> None:
    """Leakage-capped recompression pressure vs the no-leak bound (PHYSICS.md 3.7)."""

    geometry = curve.geometry
    cycle = leaky_cycle(geometry, trace=curve.trace)
    leaky_deg = [degrees(a) for a in cycle.angles_rad]
    leaky_mpa = [p * PA_TO_MPA for p in cycle.compression_pressure_pa]

    fig, ax = plt.subplots(figsize=(8.4, 5.4), dpi=140)
    ax.plot(
        curve.angles_deg,
        curve.compression_mpa,
        color="#9aa0a8",
        linewidth=1.6,
        linestyle="--",
        label="누설 무시 (§3.4, 상한)",
    )
    ax.plot(leaky_deg, leaky_mpa, color=COMPRESSION_COLOR, linewidth=2.2, label="누설 반영 (§3.7)")
    ax.axhline(curve.discharge_port_mpa, color="#9aa0a8", linestyle=":", linewidth=1.0)
    ax.set_xlabel("크랭크 각 θ (deg)")
    ax.set_ylabel("압축 챔버 압력 p (MPa)")
    ax.set_title("재압축 스파이크 상한 — 누설이 완만하게 억제 (§3.7)", fontsize=10, color="none")
    ax.set_xlim(0, 360)
    ax.set_xticks(range(0, 361, 90))
    ax.grid(True, color="#e2e5ea", linewidth=0.7)
    _legend_above(ax, ncol=2)
    ax.text(
        0.97,
        0.60,
        f"누설 무시 피크 {cycle.no_leak_peak_pa * PA_TO_MPA:.1f} MPa\n"
        f"누설 반영 피크 {cycle.capped_peak_pa * PA_TO_MPA:.1f} MPa\n"
        "(재압축 구간 ~7° 로 짧아 억제 완만;\n"
        " 실제 상한은 리드밸브 역류 — 별도 모델)",
        transform=ax.transAxes,
        va="top",
        ha="right",
        fontsize=8.5,
        bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "#c7ccd4"},
    )
    _save(fig, path)


def render_chamber_mass(curve: PressureCurve, path: Path) -> None:
    """Compression-chamber mass history and volumetric efficiency (PHYSICS.md 3.7)."""

    geometry = curve.geometry
    cycle = leaky_cycle(geometry, trace=curve.trace)
    clearance_only = leaky_cycle(geometry, gap_m=1.0e-9, trace=curve.trace)
    leaky_deg = [degrees(a) for a in cycle.angles_rad]
    mass_mg = [m * 1.0e6 for m in cycle.mass_kg]

    fig, ax = plt.subplots(figsize=(8.4, 5.4), dpi=140)
    ax.plot(leaky_deg, mass_mg, color=SUCTION_COLOR, linewidth=2.2)
    ax.axhline(0.0, color="#c7ccd4", linewidth=0.8)
    ax.set_xlabel("크랭크 각 θ (deg)")
    ax.set_ylabel("챔버 질량 m (mg)")
    ax.set_title("압축 챔버 질량 이력 — 토출 배출 + 누설 (§3.7)", fontsize=10, color="none")
    ax.set_xlim(leaky_deg[0], 360)
    ax.grid(True, color="#e2e5ea", linewidth=0.7)
    ax.text(
        0.03,
        0.05,
        f"체적효율 η_v = {cycle.volumetric_efficiency:.3f}\n"
        f"  (사체적 재팽창만: {clearance_only.volumetric_efficiency:.3f})\n"
        f"토출 질량 {cycle.delivered_mass_kg * 1.0e6:.0f} mg/rev, "
        f"누설 {cycle.leaked_mass_kg * 1.0e6:.1f} mg\n"
        "가정: 유효 틈 5 μm·Cd 0.6, R410A ρ 31.4·γ 1.10 (CoolProp)",
        transform=ax.transAxes,
        va="bottom",
        ha="left",
        fontsize=8.5,
        bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "#c7ccd4"},
    )
    _save(fig, path)


# --------------------------------------------------------------------------
# Figure: discharge reed valve
# --------------------------------------------------------------------------


def render_reed_valve_pressure(curve: PressureCurve, path: Path) -> None:
    """Discharge overpressure from the finite reed-valve port area (PHYSICS.md 3.8)."""

    geometry = curve.geometry
    trace = curve.trace
    # A fine crank-angle step keeps the near-port-close forward-Euler trace smooth
    # (the peak and loss are already converged at the default resolution).
    cycle = valved_cycle(geometry, trace=trace, samples=2880)
    valved_deg = [degrees(a) for a in cycle.angles_rad]
    valved_mpa = [p * PA_TO_MPA for p in cycle.compression_pressure_pa]

    fig, ax = plt.subplots(figsize=(10.2, 6.0), dpi=140)
    ax.plot(
        curve.angles_deg,
        curve.compression_mpa,
        color="#9aa0a8",
        linewidth=1.6,
        linestyle="--",
        label="밸브 무시 (§3.4, p_dis 클램프)",
    )
    ax.plot(valved_deg, valved_mpa, color=COMPRESSION_COLOR, linewidth=2.2, label="리드밸브 (§3.8)")
    ax.axhline(curve.discharge_port_mpa, color="#9aa0a8", linestyle=":", linewidth=1.0)
    ax.set_xlabel("크랭크 각 θ (deg)")
    ax.set_ylabel("압축 챔버 압력 p (MPa)")
    ax.set_title("토출 리드밸브 — 유한 포트 면적이 만드는 과압 (§3.8)", fontsize=10, color="none")
    ax.set_xlim(0, 360)
    ax.set_xticks(range(0, 361, 90))
    ax.grid(True, color="#e2e5ea", linewidth=0.7)
    _legend_above(ax, ncol=2)
    ax.text(
        0.03,
        0.97,
        f"토출 과압 피크 {cycle.delivery_peak_pa * PA_TO_MPA:.2f} MPa "
        f"(p_dis {curve.discharge_port_mpa:.2f})\n"
        f"밸브 과압 손실 {cycle.overpressure_power_w:.0f} W (지시의 ~7%)\n"
        f"재압축 피크 {cycle.recompression_peak_pa * PA_TO_MPA:.1f} MPa - 포트 닫힘 후\n"
        " 면적 0 이라 밸브가 못 억제 (누설만, §3.7)\n"
        "준정적 -> 역류 없음 (역류는 리드 동역학 = Level 3)",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8.5,
        bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "#c7ccd4"},
    )
    _save(fig, path)


def render_reed_valve_port_area(curve: PressureCurve, path: Path) -> None:
    """Geometric discharge-port open area vs crank angle (PHYSICS.md 3.8)."""

    geometry = curve.geometry
    degrees_axis = list(range(0, 361))
    area_mm2 = [port_open_area_m2(geometry, radians(d)) * 1.0e6 for d in degrees_axis]
    window = discharge_window(geometry)

    fig, ax = plt.subplots(figsize=(10.2, 5.4), dpi=140)
    ax.plot(degrees_axis, area_mm2, color=SUCTION_COLOR, linewidth=2.2, label="토출 포트 개방 면적")
    ax.axvspan(
        degrees(window.start_rad),
        degrees(window.end_rad),
        color=MERGED_COLOR,
        alpha=0.22,
        label=f"토출 창 [{degrees(window.start_rad):.0f}, {degrees(window.end_rad):.0f}]°",
    )
    ax.set_xlabel("크랭크 각 θ (deg)")
    ax.set_ylabel("포트 개방 면적 (mm²)")
    ax.set_title("토출 포트 개방 면적 — 접촉점이 덮으며 0으로 (§3.8)", fontsize=10, color="none")
    ax.set_xlim(0, 360)
    ax.set_xticks(range(0, 361, 90))
    ax.grid(True, color="#e2e5ea", linewidth=0.7)
    _legend_above(ax, ncol=2)
    _save(fig, path)


# --------------------------------------------------------------------------
# Figure: geometric angle definitions (crank angle theta, rotor swing phi)
# --------------------------------------------------------------------------


def render_crank_angle_definition(geometry: RotaryGeometry, path: Path) -> None:
    """Cross-section defining the crank angle theta and rotor swing phi (PHYSICS.md 4.2).

    theta is the drive (crank) angle — the bore angle of the colinear O, O_r, C
    line, measured from +y (the vane) clockwise. phi is the rotor swing (a
    dependent angle) — the tilt of the rotor mouth axis from +y, oscillating
    +/-10.37 deg with theta.
    """

    theta_deg = 120.0
    angle_rad = radians(theta_deg)
    bore_mm = geometry.cylinder_radius_m / MM
    limit = bore_mm + 13.0
    state = prescribed_state(geometry, angle_rad)
    rotor_x, rotor_y = state.rotor_center_m[0] / MM, state.rotor_center_m[1] / MM
    groove_x_m, groove_y_m = state.cutout_center_m
    groove_x, groove_y = groove_x_m / MM, groove_y_m / MM
    contact = (bore_mm * sin(angle_rad), bore_mm * cos(angle_rad))

    fig, ax = plt.subplots(figsize=(8.8, 8.8), dpi=140)
    ax.set_aspect("equal")
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.axis("off")
    fig.suptitle(
        "각도 정의 — 크랭크 각 θ와 로터 스윙 " + r"$\phi$" + " (§4.2)", fontsize=13, color="none"
    )

    # Chambers (IN, +x side of the contact; OUT, the rest).
    ax.add_patch(
        Wedge((0, 0), bore_mm, 90.0 - theta_deg, 90.0, facecolor=SUCTION_COLOR, alpha=0.16)
    )
    ax.add_patch(
        Wedge((0, 0), bore_mm, -270.0, 90.0 - theta_deg, facecolor=COMPRESSION_COLOR, alpha=0.16)
    )
    # Rotor, bush, vane, bore.
    contour = rotor_contour(geometry, angle_rad)
    ax.fill(
        [p[0] / MM for p in contour.material],
        [p[1] / MM for p in contour.material],
        facecolor=ROTOR_COLOR,
        edgecolor="black",
        linewidth=1.5,
        zorder=2,
    )
    bush = SwingBush()
    for side in (1.0, -1.0):
        xs, ys = _bush_outline_mm(groove_x_m + side * bush.piece_shift_m, groove_y_m, side)
        ax.fill(xs, ys, facecolor=BUSH_COLOR, edgecolor="black", linewidth=0.8, zorder=3)
    vane_x, vane_y = _vane_outline_mm(geometry, state.vane_tip_m[1])
    ax.fill(vane_x, vane_y, facecolor=VANE_COLOR, edgecolor="black", linewidth=1.2, zorder=4)
    ax.add_patch(
        plt.Circle((0, 0), bore_mm, facecolor="none", edgecolor="black", linewidth=2.0, zorder=5)
    )

    # Coordinate frame at O: x toward inlet, y toward the vane.
    ax.annotate(
        "",
        xy=(bore_mm * 0.55, 0),
        xytext=(0, 0),
        arrowprops={"arrowstyle": "-|>", "color": "#333", "lw": 1.8},
        zorder=7,
    )
    ax.text(bore_mm * 0.57, -2.6, "x (인렛)", fontsize=12, color="#333", zorder=8)
    ax.annotate(
        "",
        xy=(0, 11.5),
        xytext=(0, 0),
        arrowprops={"arrowstyle": "-|>", "color": "#333", "lw": 1.8},
        zorder=7,
    )
    ax.text(-1.5, 12.2, "y (베인)", fontsize=12, color="#333", ha="right", zorder=8)

    # Colinear O -> O_r -> C line, the crank direction at angle theta.
    ax.plot([0, contact[0]], [0, contact[1]], color="#d62728", linewidth=1.3, zorder=6)
    ax.plot(0, 0, "+", color="black", markersize=12, markeredgewidth=2, zorder=8)
    ax.text(-2.5, -3.0, "O (크랭크 축)", fontsize=11, ha="right", zorder=8)
    ax.plot(rotor_x, rotor_y, "+", color=GUIDE_COLOR, markersize=11, markeredgewidth=2, zorder=8)
    ax.text(rotor_x + 1.6, rotor_y - 2.2, "O_r (로터 중심)", fontsize=11, zorder=8)
    ax.plot(*contact, "o", color="#d62728", markersize=7, zorder=8)
    ax.text(
        contact[0] - 1.0,
        contact[1] - 2.2,
        "C (로터-실린더 접촉)",
        fontsize=11,
        ha="right",
        zorder=8,
    )

    # theta arc at O: from +y to the O->C direction, clockwise.
    arc_r = bore_mm * 0.42
    ax.add_patch(
        Arc(
            (0, 0),
            2 * arc_r,
            2 * arc_r,
            angle=0.0,
            theta1=90.0 - theta_deg,
            theta2=90.0,
            color="#7a4fb0",
            linewidth=2.4,
            zorder=7,
        )
    )
    mid = radians(90.0 - theta_deg * 0.5)
    ax.text(
        arc_r * 1.12 * cos(mid),
        arc_r * 1.12 * sin(mid),
        "θ",
        color="#7a4fb0",
        fontsize=17,
        ha="center",
        va="center",
        zorder=8,
    )

    # phi (rotor swing) at O_r: +y reference vs the mouth axis O_r -> groove.
    ax.plot(
        [rotor_x, rotor_x],
        [rotor_y, rotor_y + 12.0],
        color="#2e7d32",
        linestyle=":",
        linewidth=1.2,
        zorder=6,
    )
    ax.plot([rotor_x, groove_x], [rotor_y, groove_y], color="#2e7d32", linewidth=1.5, zorder=6)
    swing = degrees(atan2(groove_x - rotor_x, groove_y - rotor_y))
    parc_r = 9.0
    lo, hi = sorted((90.0, 90.0 - swing))
    ax.add_patch(
        Arc(
            (rotor_x, rotor_y),
            2 * parc_r,
            2 * parc_r,
            angle=0.0,
            theta1=lo,
            theta2=hi,
            color="#2e7d32",
            linewidth=2.2,
            zorder=7,
        )
    )
    ax.text(
        rotor_x - 1.5,
        rotor_y + parc_r + 1.6,
        r"$\phi$",
        color="#2e7d32",
        fontsize=17,
        ha="center",
        zorder=8,
    )

    ax.text(
        0.02,
        0.02,
        "θ: 크랭크(구동) 각 — O·O_r·접촉점 C 일직선의 방향, +y(베인)에서 시계방향\n"
        r"$\phi$: 로터 스윙(종속) — 마우스 축이 베인(+y)에서 기운 각, θ에 따라 ±10.37°"
        "\ne: O–O_r 편심 4.5 mm · IN(흡입, +x쪽)/OUT(압축·토출)는 접촉점 C에서 분리",
        transform=ax.transAxes,
        va="bottom",
        ha="left",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "#c7ccd4"},
        zorder=9,
    )
    _save(fig, path)


# --------------------------------------------------------------------------
# Figure: gas force / bearing reaction free-body with axes and torque axis
# --------------------------------------------------------------------------


def render_force_diagram(curve: PressureCurve, path: Path) -> None:
    """Compressor free body with the coordinate + torque axes (PHYSICS.md 4.2/4.5/4.6)."""

    geometry = curve.geometry
    theta_deg = 235.0
    angle_rad = radians(theta_deg)
    bore_mm = geometry.cylinder_radius_m / MM
    limit = bore_mm + 13.0
    state = prescribed_state(geometry, angle_rad)

    fig, ax = plt.subplots(figsize=(8.6, 8.6), dpi=140)
    ax.set_aspect("equal")
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.axis("off")

    # Suction (IN, +x) and discharge (OUT, -x) chambers, split at the contact.
    ax.add_patch(
        Wedge((0.0, 0.0), bore_mm, 90.0 - theta_deg, 90.0, facecolor=SUCTION_COLOR, alpha=0.22)
    )
    ax.add_patch(
        Wedge(
            (0.0, 0.0), bore_mm, -270.0, 90.0 - theta_deg, facecolor=COMPRESSION_COLOR, alpha=0.22
        )
    )

    contour = rotor_contour(geometry, angle_rad)
    ax.fill(
        [p[0] / MM for p in contour.material],
        [p[1] / MM for p in contour.material],
        facecolor=ROTOR_COLOR,
        edgecolor="black",
        linewidth=1.5,
        zorder=2,
    )
    groove_x, groove_y = state.cutout_center_m
    bush = SwingBush()
    for side in (1.0, -1.0):
        xs, ys = _bush_outline_mm(groove_x + side * bush.piece_shift_m, groove_y, side)
        ax.fill(xs, ys, facecolor=BUSH_COLOR, edgecolor="black", linewidth=0.8, zorder=3)
    vane_x, vane_y = _vane_outline_mm(geometry, state.vane_tip_m[1])
    ax.fill(vane_x, vane_y, facecolor=VANE_COLOR, edgecolor="black", linewidth=1.2, zorder=4)
    ax.add_patch(
        plt.Circle(
            (0.0, 0.0), bore_mm, facecolor="none", edgecolor="black", linewidth=2.0, zorder=5
        )
    )
    ax.plot(
        bore_mm * sin(angle_rad),
        bore_mm * cos(angle_rad),
        "o",
        color="#c0392b",
        markersize=7,
        zorder=6,
    )

    rotor_x, rotor_y = state.rotor_center_m[0] / MM, state.rotor_center_m[1] / MM

    # Coordinate frame at O: x toward inlet (+x), y toward the vane (+y).
    axis_len = bore_mm * 0.52
    for (tx, ty), label, (lx, ly) in (
        ((axis_len, 0.0), "x", (axis_len + 2.5, -2.5)),
        ((0.0, axis_len), "y", (2.0, axis_len + 2.0)),
    ):
        ax.annotate(
            "",
            xy=(tx, ty),
            xytext=(0.0, 0.0),
            arrowprops={"arrowstyle": "-|>", "color": "black", "lw": 2.2},
            zorder=8,
        )
        ax.text(lx, ly, label, fontsize=15, ha="center", va="center", zorder=9)
    # Torque axis z out of plane at O.
    ax.add_patch(
        plt.Circle((0.0, 0.0), 2.4, facecolor="white", edgecolor="black", linewidth=1.6, zorder=9)
    )
    ax.plot(0.0, 0.0, "o", color="black", markersize=3.5, zorder=10)
    ax.text(3.0, 2.8, "z", fontsize=15, ha="center", va="center", zorder=9)
    # Torque sense about z (curved moment arrow, discharge side) labelled T_gas.
    t_r = 12.0
    ax.annotate(
        "",
        xy=(t_r * cos(radians(235)), t_r * sin(radians(235))),
        xytext=(t_r * cos(radians(125)), t_r * sin(radians(125))),
        arrowprops={
            "arrowstyle": "-|>",
            "color": "#7a4fb0",
            "lw": 2.4,
            "connectionstyle": "arc3,rad=0.32",
        },
        zorder=8,
    )
    ax.text(-14.8, -3.5, "T_gas", fontsize=14, color="#7a4fb0", ha="center", zorder=9)

    # Gas force on the rotor (mouth-aware, §4.5) and the equal-opposite reaction.
    load = true_gas_load(geometry, angle_rad, trace=curve.trace)
    fx, fy = load.rotor_force_n
    force_scale = (bore_mm * 0.42) / max(hypot(fx, fy), 1.0)
    ax.plot(rotor_x, rotor_y, "+", color=GUIDE_COLOR, markersize=10, markeredgewidth=2, zorder=10)
    ax.annotate(
        "",
        xy=(rotor_x + fx * force_scale, rotor_y + fy * force_scale),
        xytext=(rotor_x, rotor_y),
        arrowprops={"arrowstyle": "-|>", "color": COMPRESSION_COLOR, "lw": 3.0},
        zorder=11,
    )
    ax.text(
        rotor_x + fx * force_scale + 1.5,
        rotor_y + fy * force_scale - 1.5,
        "F_gas",
        fontsize=14,
        color=COMPRESSION_COLOR,
        zorder=12,
    )
    ax.annotate(
        "",
        xy=(rotor_x - fx * force_scale, rotor_y - fy * force_scale),
        xytext=(rotor_x, rotor_y),
        arrowprops={"arrowstyle": "-|>", "color": SUCTION_COLOR, "lw": 3.0},
        zorder=11,
    )
    ax.text(
        rotor_x - fx * force_scale - 1.5,
        rotor_y - fy * force_scale + 1.5,
        "R_j",
        fontsize=14,
        color=SUCTION_COLOR,
        ha="right",
        zorder=12,
    )
    _reed_note(ax, x=0.02, y=0.02, va="bottom")
    _save(fig, path)


# --------------------------------------------------------------------------
# Figure: indicated-work thermodynamic cross-check
# --------------------------------------------------------------------------


def render_thermo_crosscheck(curve: PressureCurve, path: Path) -> None:
    """Indicated work: P-V loop vs CoolProp isentropic route (PHYSICS.md 3.5)."""

    geometry = curve.geometry
    trace = curve.trace
    check = isentropic_cross_check(geometry, trace=trace)
    displacement_m3 = 22.02e-6  # total gas at TDC
    ideal_full_w = (
        SUCTION_DENSITY_KG_M3
        * displacement_m3
        * check.isentropic_enthalpy_rise_j_kg
        * check.frequency_hz
    )

    labels = [
        "P-V 적분 (경로 A, §3.5)",
        "열역학 m·Δh_s (경로 B, CoolProp)",
        "등엔트로피 전변위 (상한, η_v=1)",
    ]
    values = [check.indicated_power_w, check.isentropic_power_w, ideal_full_w]
    colors = [COMPRESSION_COLOR, SUCTION_COLOR, "#c7ccd4"]

    fig, ax = plt.subplots(figsize=(10.2, 4.8), dpi=140)
    positions = range(len(values))
    ax.barh(list(positions), values, color=colors, height=0.62)
    ax.set_yticks(list(positions))
    ax.set_yticklabels(labels, fontsize=9.5)
    ax.invert_yaxis()
    for position, value in zip(positions, values, strict=True):
        ax.text(value + 6.0, position, f"{value:.0f} W", va="center", fontsize=10)
    ax.set_xlabel("지시 / 등가 동력 (W)")
    ax.set_xlim(0, max(values) * 1.18)
    ax.grid(True, axis="x", color="#e2e5ea", linewidth=0.7)
    ax.set_axisbelow(True)
    # Title kept in the file but invisible (filename is the identity, PPT adds its own).
    ax.set_title("지시 일 열역학 교차검증 — P-V vs 등엔트로피 (§3.5)", fontsize=10, color="none")
    ax.text(
        0.98,
        0.06,
        f"비지시일 W_ind/m = {check.specific_work_j_kg / 1e3:.1f} kJ/kg\n"
        f"  ~ CoolProp Δh_s = {check.isentropic_enthalpy_rise_j_kg / 1e3:.2f} kJ/kg "
        f"(오차 {check.relative_error * 100:.1f}%)\n"
        "독립 2경로 일치 -> 실측 P-V 없이 지시 일 검증\n"
        "(등엔트로피-가역 하한; 실측 인디케이터는 +15~40%)",
        transform=ax.transAxes,
        va="bottom",
        ha="right",
        fontsize=8.5,
        bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "#c7ccd4"},
    )
    _save(fig, path)


# --------------------------------------------------------------------------
# Pruning and entry point
# --------------------------------------------------------------------------


# Registry of every figure this script owns. ``kind`` is "curve" when the
# renderer needs the shared port-timed pressure trace and "geometry" when the
# static geometry is enough. Add a renderer here (and it is regenerated and
# kept by --prune); legacy figures are ported into this table one at a time.
def _reynolds_curved_data(geometry: RotaryGeometry) -> dict:
    """Bush curved film: 1-D ``arc_film`` load vs the long-bearing reference over ε.

    Shared by the ``reynolds_curved`` figure and its ``.dat`` export.
    """

    bush = SwingBush()
    radius = bush.piece_outer_radius_m
    height = geometry.cylinder_height_m - 2.0 * AXIAL_CLEARANCE_M
    curved_gap = geometry.cutout_radius_m - radius
    entrain, visc = 50.0, LUBRICANT_VISCOSITY_PA_S
    eps = [0.1 + 0.05 * i for i in range(17)]
    arc_mag, lb_mag = [], []
    for e in eps:
        arc = arc_film_force(
            -e * curved_gap,
            0.0,
            0.0,
            0.0,
            entrain,
            arc_center_rad=pi / 2,
            arc_half_span_rad=pi / 2,
            radius_m=radius,
            length_m=height,
            clearance_m=curved_gap,
            viscosity_pa_s=visc,
            n_beta=2001,
        )
        ref = long_bearing_load(
            e,
            entrain,
            radius_m=radius,
            length_m=height,
            clearance_m=curved_gap,
            viscosity_pa_s=visc,
            condition="half",
        )
        arc_mag.append(hypot(arc.force_x_n, arc.force_y_n))
        lb_mag.append(ref.magnitude_n)
    err = max(abs(a - b) / b for a, b in zip(arc_mag, lb_mag, strict=True))
    return {"eps": eps, "arc_mag": arc_mag, "lb_mag": lb_mag, "err": err}


def _reynolds_flat_data(geometry: RotaryGeometry) -> dict:
    """Bush flat film: 1-D ``slider_film`` load vs the fixed-incline slider over the tilt."""

    height = geometry.cylinder_height_m - 2.0 * AXIAL_CLEARANCE_M
    visc = LUBRICANT_VISCOSITY_PA_S
    lf, cf, u = 11.0e-3, 10.0e-6, 0.86
    tilts = [1.0e-4 + (5.0e-4) * i / 15 for i in range(16)]
    sld, inc = [], []
    for g in tilts:
        f = flat_slider_film(
            0.0, -g, 0.0, 0.0, u, length_m=lf, height_m=height, clearance_m=cf, n_s=4001
        )
        h1, h0 = cf + 0.5 * lf * g, cf - 0.5 * lf * g
        a = h1 / h0
        w = (
            (6.0 * visc * u * lf**2)
            / ((a - 1.0) ** 2 * h0**2)
            * (log(a) - 2.0 * (a - 1.0) / (a + 1.0))
        )
        sld.append(f.normal_force_n)
        inc.append(w * height)
    err = max(abs(a - b) / b for a, b in zip(sld, inc, strict=True))
    return {"tilts": tilts, "sld": sld, "inc": inc, "err": err}


def _reynolds_journal_data() -> dict:
    """Journal film: 1-D numerical Reynolds load vs the Ocvirk short bearing over ε."""

    epsj = [0.05 + 0.05 * i for i in range(18)]
    ocv = [short_bearing_force(e, 95.0).magnitude_n for e in epsj]
    num = [solve_short_bearing_1d(e, 95.0).magnitude_n for e in epsj]
    err = max(abs(n - o) / o for n, o in zip(num, ocv, strict=True))
    return {"epsj": epsj, "ocv": ocv, "num": num, "err": err}


def render_film_reynolds_curved(geometry: RotaryGeometry, path: Path) -> None:
    """Bush curved film: 1-D ``arc_film`` vs the long-bearing Sommerfeld load (§4.11).

    The axial-uniform arc-film reduction reproduces the closed-form long-bearing load to
    ~1e-4 across the eccentricity range. (One of three 1-D-Reynolds validity figures — see
    also ``reynolds_flat_vs_incline_slider`` and ``reynolds_journal_vs_ocvirk``.)
    """

    blue, red = "#2166ac", "#c0392b"
    d = _reynolds_curved_data(geometry)
    eps, arc_mag, lb_mag, err_curved = d["eps"], d["arc_mag"], d["lb_mag"], d["err"]

    fig, ax = plt.subplots(figsize=(8.4, 5.6), dpi=140)
    ax.grid(True, which="both", color="#e2e5ea", lw=0.7)
    ax.semilogy(eps, lb_mag, "-", color=blue, lw=2.2, label="long-bearing 해석해 (§4.11)")
    ax.semilogy(
        eps, arc_mag, "o", color=red, ms=6, mfc="none", mew=1.6, label="1-D Reynolds arc_film"
    )
    ax.set_xlabel("편심비 ε")
    ax.set_ylabel("|F| (N)")
    ax.set_title(
        f"부시 곡면 유막 — arc_film vs long-bearing (§4.11)\n최대 상대오차 {err_curved:.1e}",
        fontsize=12,
        fontweight="bold",
    )
    ax.legend(fontsize=9)
    _save(fig, path)


def render_film_reynolds_flat(geometry: RotaryGeometry, path: Path) -> None:
    """Bush flat film: 1-D ``slider_film`` vs the Reynolds fixed-incline slider (§4.11).

    The 1-D pad reduction reproduces the closed-form fixed-incline slider load to ~1e-4
    across the wedge-tilt range. (See also ``reynolds_curved_vs_long_bearing`` and
    ``reynolds_journal_vs_ocvirk``.)
    """

    blue, red = "#2166ac", "#c0392b"
    d = _reynolds_flat_data(geometry)
    tilts, sld, inc, err_flat = d["tilts"], d["sld"], d["inc"], d["err"]

    fig, ax = plt.subplots(figsize=(8.4, 5.6), dpi=140)
    ax.grid(True, which="both", color="#e2e5ea", lw=0.7)
    ax.plot([t * 1e3 for t in tilts], inc, "-", color=blue, lw=2.2, label="Reynolds 경사판 해석해")
    ax.plot(
        [t * 1e3 for t in tilts],
        sld,
        "o",
        color=red,
        ms=6,
        mfc="none",
        mew=1.6,
        label="1-D Reynolds slider_film",
    )
    ax.set_xlabel("쐐기 기울기 |γ| (×1e-3 rad)")
    ax.set_ylabel("법선 하중 (N)")
    ax.set_title(
        f"부시 평면 유막 — slider_film vs 경사판 (§4.11)\n최대 상대오차 {err_flat:.1e}",
        fontsize=12,
        fontweight="bold",
    )
    ax.legend(fontsize=9)
    _save(fig, path)


def render_film_reynolds_journal(geometry: RotaryGeometry, path: Path) -> None:
    """Journal film: 1-D numerical Reynolds vs the Ocvirk short bearing (§4.12).

    The 1-D short-bearing solve reproduces the closed-form Ocvirk load to ~1e-4 across the
    eccentricity range. (See also ``reynolds_curved_vs_long_bearing`` and
    ``reynolds_flat_vs_incline_slider``.)
    """

    blue, red = "#2166ac", "#c0392b"
    d = _reynolds_journal_data()
    epsj, ocv, num, err_jrnl = d["epsj"], d["ocv"], d["num"], d["err"]

    fig, ax = plt.subplots(figsize=(8.4, 5.6), dpi=140)
    ax.grid(True, which="both", color="#e2e5ea", lw=0.7)
    ax.semilogy(epsj, ocv, "-", color=blue, lw=2.2, label="Ocvirk 단축 해석해 (§4.9)")
    ax.semilogy(
        epsj, num, "o", color=red, ms=6, mfc="none", mew=1.6, label="1-D 수치 Reynolds (§4.12)"
    )
    ax.set_xlabel("편심비 ε")
    ax.set_ylabel("|F| (N)")
    ax.set_title(
        f"저널 유막 — 1-D Reynolds vs Ocvirk (§4.12)\n최대 상대오차 {err_jrnl:.1e}",
        fontsize=12,
        fontweight="bold",
    )
    ax.legend(fontsize=9)
    _save(fig, path)


def _coupled_bush_orbit(curve: PressureCurve):
    """The coupled 9-DOF rotor+bush orbit (seal_contact=True), computed once and memoised on
    the curve so the clearance/friction figures share one (expensive) integration.

    Stored on the ``curve`` instance rather than a module-global ``id()`` map: the cache then
    lives exactly as long as the curve, with no risk of a freed curve's id being reused for a
    different geometry.
    """

    orbit = getattr(curve, "_coupled_orbit", None)
    if orbit is None:
        print("  integrating the coupled 9-DOF rotor+bush orbit (a few minutes)...")
        orbit = integrate_rotor_bush_orbit(
            curve.geometry,
            revolutions=4,
            samples=180,
            grid_samples=180,
            n_beta=121,
            n_s=81,
            trace=curve.trace,
        )
        curve._coupled_orbit = orbit
    return orbit


def _clearance_at_angle_state(curve: PressureCurve) -> dict:
    """Coupled-orbit geometry at the crank angle where the bush curved film is thinnest.

    Shared by the three ``film_clearance_*`` figures (journal, bush curved, bush flat) so the
    most-loaded instant and its piece/groove positions are computed once per curve.
    """

    geometry = curve.geometry
    orbit = _coupled_bush_orbit(curve)
    bush = SwingBush()
    flat_gap, _ = film_thicknesses_m(geometry, bush)
    lever = geometry.cutout_offset_m
    throw = geometry.eccentricity_m

    jmin = min(
        range(len(orbit.crank_angle_rad)),
        key=lambda i: min(orbit.in_curved_film_m[i], orbit.out_curved_film_m[i]),
    )
    theta = orbit.crank_angle_rad[jmin]
    phi_orient = (
        prescribed_state(geometry, theta).rotor_orientation_rad
        + orbit.rotor_attitude_deviation_rad[jmin]
    )
    e_jx, e_jy = orbit.eccentricity_x_m[jmin], orbit.eccentricity_y_m[jmin]
    pin = (throw * sin(theta), throw * cos(theta))
    bore = (pin[0] + e_jx, pin[1] + e_jy)
    groove = (bore[0] + lever * cos(phi_orient), bore[1] + lever * sin(phi_orient))
    pieces = (
        (
            1.0,
            orbit.in_piece_x_m[jmin],
            orbit.in_piece_y_m[jmin],
            orbit.in_piece_attitude_rad[jmin],
            "IN",
        ),
        (
            -1.0,
            orbit.out_piece_x_m[jmin],
            orbit.out_piece_y_m[jmin],
            orbit.out_piece_attitude_rad[jmin],
            "OUT",
        ),
    )
    return {
        "geometry": geometry,
        "bush": bush,
        "half_arc": bush.half_arc_rad(),
        "curved_gap": geometry.cutout_radius_m - bush.piece_outer_radius_m,
        "flat_gap": flat_gap,
        "shift": bush.piece_shift_m,
        "c_j": orbit.journal_clearance_m,
        "theta": theta,
        "e_jx": e_jx,
        "e_jy": e_jy,
        "groove": groove,
        "lf": flat_contact_length_m(geometry, theta, bush),
        "pieces": pieces,
    }


def render_film_clearance_journal(curve: PressureCurve, path: Path) -> None:
    """Journal film thickness h(φ) around the circumference at the most-loaded angle (§4.14)."""

    st = _clearance_at_angle_state(curve)
    c_j, e_jx, e_jy = st["c_j"], st["e_jx"], st["e_jy"]
    blue, red, grey = "#2166ac", "#c0392b", "#8a8f98"

    fig, ax = plt.subplots(figsize=(8.4, 5.6), dpi=140)
    ax.grid(True, color="#e2e5ea", lw=0.7)
    phi = [2.0 * pi * i / 360 for i in range(361)]
    h_j = [c_j - (e_jx * cos(p) + e_jy * sin(p)) for p in phi]
    j_arg = min(range(len(h_j)), key=lambda i: h_j[i])
    ax.plot([degrees(p) for p in phi], [h * 1e6 for h in h_j], "-", color=blue, lw=2.2)
    ax.axhline(c_j * 1e6, ls="--", color=grey, lw=1.3, label=f"공칭 c_j={c_j * 1e6:.0f} µm")
    ax.plot(degrees(phi[j_arg]), h_j[j_arg] * 1e6, "v", color=red, ms=9)
    ax.set_xlabel("원주각 φ (deg)")
    ax.set_ylabel("유막 두께 h (µm)")
    ax.set_xlim(0, 360)
    ax.set_title(
        f"저널 클리어런스 (θ={degrees(st['theta']):.0f}°) — "
        f"min {min(h_j) * 1e6:.1f} µm (ε={hypot(e_jx, e_jy) / c_j:.2f})",
        fontsize=12,
        fontweight="bold",
    )
    ax.legend(fontsize=9)
    _save(fig, path)


def render_film_clearance_bush_curved(curve: PressureCurve, path: Path) -> None:
    """Bush curved (rotor-groove) film h(β) for both pieces at the most-loaded angle (§4.14)."""

    st = _clearance_at_angle_state(curve)
    half_arc, curved_gap, groove = st["half_arc"], st["curved_gap"], st["groove"]
    grey = "#8a8f98"

    fig, ax = plt.subplots(figsize=(8.4, 5.6), dpi=140)
    ax.grid(True, color="#e2e5ea", lw=0.7)
    for side, x_k, y_k, phi_k, tag in st["pieces"]:
        ecc_x, ecc_y = x_k - groove[0], y_k - groove[1]
        center = (0.0 if side > 0 else pi) + phi_k
        beta = [center - half_arc + 2.0 * half_arc * i / 360 for i in range(361)]
        h_c = [curved_gap - (ecc_x * cos(b) + ecc_y * sin(b)) for b in beta]
        ax.plot(
            [degrees(b - center) for b in beta],
            [h * 1e6 for h in h_c],
            lw=2.2,
            label=f"{tag} (min {min(h_c) * 1e6:.2f} µm)",
        )
    ax.axhline(
        curved_gap * 1e6, ls="--", color=grey, lw=1.3, label=f"공칭 gap={curved_gap * 1e6:.0f} µm"
    )
    ax.set_xlabel("아크각 (β - 중심) (deg)")
    ax.set_ylabel("유막 두께 h (µm)")
    ax.set_title(
        f"부시 곡면 클리어런스 (두 조각, θ={degrees(st['theta']):.0f}°)",
        fontsize=12,
        fontweight="bold",
    )
    ax.legend(fontsize=9)
    _save(fig, path)


def render_film_clearance_bush_flat(curve: PressureCurve, path: Path) -> None:
    """Bush flat (vane) film thickness h(s) for both pieces at the loaded angle (§4.14)."""

    st = _clearance_at_angle_state(curve)
    flat_gap, shift, lf = st["flat_gap"], st["shift"], st["lf"]
    grey = "#8a8f98"

    fig, ax = plt.subplots(figsize=(8.4, 5.6), dpi=140)
    ax.grid(True, color="#e2e5ea", lw=0.7)
    s_grid = [-0.5 * lf + lf * i / 200 for i in range(201)]
    for side, x_k, _y_k, phi_k, tag in st["pieces"]:
        approach = shift - side * x_k
        h_f = [(flat_gap - approach) + s * phi_k for s in s_grid]
        ax.plot(
            [s * 1e3 for s in s_grid],
            [h * 1e6 for h in h_f],
            lw=2.2,
            label=f"{tag} (min {min(h_f) * 1e6:.2f} µm)",
        )
    ax.axhline(
        flat_gap * 1e6, ls="--", color=grey, lw=1.3, label=f"공칭 c_f={flat_gap * 1e6:.0f} µm"
    )
    ax.set_xlabel("패드 위치 s (mm)")
    ax.set_ylabel("유막 두께 h (µm)")
    ax.set_title(
        f"부시 평면 클리어런스 (두 조각, θ={degrees(st['theta']):.0f}°)",
        fontsize=12,
        fontweight="bold",
    )
    ax.legend(fontsize=9)
    _save(fig, path)


def render_friction_dynamic_vs_quasistatic(curve: PressureCurve, path: Path) -> None:
    """Per-film friction loss: quasi-static vs dynamic/coupled (§4.14/4.15).

    Bush and rotor-cylinder seal roughly triple from quasi-static to dynamic/coupled
    (thinner loaded films; vane-constraint load transfer); the journal is nearly
    unchanged (the squeeze film attenuates the peak eccentricity). The seal row is the
    coupled 9-DOF value (~11 W) vs the rigid quasi-static N_c estimate.
    """

    geometry = curve.geometry
    orbit = _coupled_bush_orbit(curve)
    rotor = integrate_rotor_orbit(geometry, trace=curve.trace)
    qs_seal = rotor_cylinder_friction_power_w(
        geometry, lambda th: contact_normal_force_n(geometry, th, trace=curve.trace)
    )
    rows = [
        (
            "부시\n(곡면+평면)",
            orbit.quasi_static_bush_friction_power_w,
            orbit.bush_friction_power_w,
        ),
        (
            "크랭크핀 저널\n(편심)",
            rotor.quasi_static_journal_friction_w,
            rotor.dynamic_journal_friction_w,
        ),
        ("로터-실린더\n실링", qs_seal, orbit.seal_contact_friction_power_w),
    ]
    fig, ax = plt.subplots(figsize=(8.6, 5.2), dpi=140)
    xs = list(range(len(rows)))
    w = 0.38
    qs = [r[1] for r in rows]
    dyn = [r[2] for r in rows]
    ax.bar([x - w / 2 for x in xs], qs, w, color="#8a8f98", label="준정적")
    ax.bar([x + w / 2 for x in xs], dyn, w, color="#c0392b", label="동역학 / 결합")
    for x, (q, d) in zip(xs, zip(qs, dyn, strict=True), strict=True):
        ax.text(x - w / 2, q, f"{q:.2f}", ha="center", va="bottom", fontsize=9)
        ax.text(x + w / 2, d, f"{d:.2f}", ha="center", va="bottom", fontsize=9)
    ax.set_xticks(xs)
    ax.set_xticklabels([r[0] for r in rows])
    ax.set_ylabel("마찰 손실 (W)")
    total_qs, total_dyn = sum(qs), sum(dyn)
    # Summary in the title (above the data), not a box over the bars — the "upper left"
    # corner is empty here since the tallest bars (seal) sit on the right.
    ax.set_title(
        "유막별 마찰 손실 — 준정적 vs 동역학/결합 (§4.14/4.15)\n"
        f"합계 준정적 {total_qs:.1f} W → 동역학 {total_dyn:.1f} W (×{total_dyn / total_qs:.2f})",
        fontsize=12,
        fontweight="bold",
    )
    ax.set_ylim(0, 1.22 * max(max(qs), max(dyn)))  # headroom for the tallest bar + label
    ax.grid(True, axis="y", color="#e2e5ea", lw=0.7)
    ax.legend(fontsize=10, loc="upper left")
    _save(fig, path)


_ASM_COLORS = {"rotor": "#d9d9d9", "bush": "#aeb8c7", "vane": "#5a5a5a", "pin": "#c9b08a"}
_ASM_K = 60  # uniform clearance-exaggeration factor shared by the bush and journal views


def _assembly_state_42(curve: PressureCurve) -> dict:
    """Geometry + coupled-orbit state at the most-loaded crank angle (~42 deg), in mm."""

    geometry = curve.geometry
    orbit = _coupled_bush_orbit(curve)
    bush = SwingBush()
    ang = [degrees(a) for a in orbit.crank_angle_rad]
    j = min(range(len(ang)), key=lambda i: abs(ang[i] - 42.0))
    thr = orbit.crank_angle_rad[j]
    state = prescribed_state(geometry, thr)
    phi = state.rotor_orientation_rad + orbit.rotor_attitude_deviation_rad[j]
    mm = 1.0 / MM
    throw = geometry.eccentricity_m * mm
    pin = np.array([throw * sin(thr), throw * cos(thr)])
    bore = pin + np.array([orbit.eccentricity_x_m[j], orbit.eccentricity_y_m[j]]) * mm
    groove = bore + geometry.cutout_offset_m * mm * np.array([cos(phi), sin(phi)])
    in_op = np.array([orbit.in_piece_x_m[j], orbit.in_piece_y_m[j]]) * mm
    out_op = np.array([orbit.out_piece_x_m[j], orbit.out_piece_y_m[j]]) * mm
    return {
        "geometry": geometry,
        "bush": bush,
        "state": state,
        "mm": mm,
        "th": ang[j],
        "Rc": geometry.cylinder_radius_m * mm,
        "Rr": geometry.rotor_radius_m * mm,
        "Rpin": 14.2,
        "cj": 0.015,
        "rcut": geometry.cutout_radius_m * mm,
        "vane_half": 0.5 * geometry.vane_width_m * mm,
        "pin": pin,
        "bore": bore,
        "groove": groove,
        "ej": np.array([orbit.eccentricity_x_m[j], orbit.eccentricity_y_m[j]]) * 1e6,
        "in_c": orbit.in_curved_film_m[j] * 1e6,
        "in_f": orbit.in_flat_film_m[j] * 1e6,
        "out_c": orbit.out_curved_film_m[j] * 1e6,
        "out_f": orbit.out_flat_film_m[j] * 1e6,
        "in_op": in_op,
        "out_op": out_op,
        "in_phi": orbit.in_piece_attitude_rad[j],
        "out_phi": orbit.out_piece_attitude_rad[j],
        "in_ecc": float(np.hypot(*(in_op - groove))) * 1e3,
        "out_ecc": float(np.hypot(*(out_op - groove))) * 1e3,
        "in_vx": in_op[0] * 1e3,
        "out_vx": out_op[0] * 1e3,
    }


def render_assembly_layout(curve: PressureCurve, path: Path) -> None:
    """Full cross-section at the most-loaded crank angle, true scale (§4.14).

    The arrangement of the journal (crank pin), rotor, fixed vane and the two swing-bush
    pieces, with the shaft axis O, crank-pin centre O_j, groove centre O_g and the seal
    contact marked. Clearances are micrometre-scale (invisible here); the bush and journal
    clearance figures show them exaggerated.
    """

    s = _assembly_state_42(curve)
    mm, Rc, Rr, Rpin, cj = s["mm"], s["Rc"], s["Rr"], s["Rpin"], s["cj"]
    rcut, pin, bore, groove = s["rcut"], s["pin"], s["bore"], s["groove"]
    c = _ASM_COLORS
    fig, ax = plt.subplots(figsize=(7.4, 7.0), dpi=140)
    ax.set_aspect("equal")
    ax.grid(True, color="#eef0f3", lw=0.6)
    ax.add_patch(plt.Circle((0, 0), Rc, fill=False, ec="#333", lw=1.6))
    ax.add_patch(plt.Circle(bore, Rr, facecolor=c["rotor"], ec="#555", lw=1.2))
    ax.add_patch(plt.Circle(bore, Rpin + cj, facecolor="white", ec="#999", lw=0.7))
    ax.add_patch(plt.Circle(pin, Rpin, facecolor=c["pin"], ec="#8a7", lw=1.0))
    vx, vy = _vane_outline_mm(s["geometry"], s["state"].vane_tip_m[1])
    ax.fill(vx, vy, facecolor=c["vane"], ec="#333", lw=1.0, zorder=5)
    for side in (1.0, -1.0):
        bx, by = _bush_outline_mm(
            groove[0] / mm + side * s["bush"].piece_shift_m, groove[1] / mm, side
        )
        ax.fill(bx, by, facecolor=c["bush"], ec="#222", lw=1.0, zorder=4)
    ax.add_patch(plt.Circle(groove, rcut, fill=False, ec="#c0392b", lw=1.0, ls=":"))
    for p, name, off in (
        (np.zeros(2), "O (축)", (1.5, -3.2)),
        (pin, "O_j", (1.8, -2.6)),
        (groove, "O_g", (2.6, 0.5)),
    ):
        ax.plot(*p, "+", color="k", ms=8, mew=1.4)
        ax.annotate(name, p, p + np.array(off), fontsize=9)
    seal = bore + Rr * bore / np.hypot(*bore)
    ax.annotate(
        "실링 접촉",
        seal,
        seal + np.array([4, 3]),
        fontsize=9,
        color="#c0392b",
        arrowprops={"arrowstyle": "->", "color": "#c0392b"},
    )
    ax.annotate("베인(고정)", (0, 34), (12, 36), fontsize=9, color="#333")
    ax.annotate(
        "두 부시",
        groove + np.array([0, 2]),
        (-20, 24),
        fontsize=9,
        color="#333",
        arrowprops={"arrowstyle": "->", "color": "#333"},
    )
    ax.set_xlim(-Rc - 4, Rc + 12)
    ax.set_ylim(-Rc - 4, Rc + 8)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title(
        f"θ={s['th']:.0f}° 전체 배치 (실척) — 저널·로터·베인·두 부시",
        fontsize=12,
        fontweight="bold",
    )
    _save(fig, path)


def render_assembly_bush_clearance(curve: PressureCurve, path: Path) -> None:
    """Bush clearances at the most-loaded angle, exaggerated x60, with the centre offsets.

    The two bush pieces flanking the fixed vane, their curved (rotor-groove) and flat
    (vane) films exaggerated by a uniform factor, plus the groove centre O_g and the two
    piece centres O_p marked relative to the vane centre-line. The IN piece is driven onto
    its groove wall (curved film into near-contact) while the OUT piece stays near-concentric.
    """

    s = _assembly_state_42(curve)
    geometry, bush, mmv = s["geometry"], s["bush"], s["mm"]
    rcut, vane_half, groove = s["rcut"], s["vane_half"], s["groove"]
    in_c, in_f, out_c, out_f = s["in_c"], s["in_f"], s["out_c"], s["out_f"]
    in_op, out_op = s["in_op"], s["out_op"]
    c, kk = _ASM_COLORS, _ASM_K
    # Rigid caricature: the piece keeps a constant shape (outer radius r_prime = rcut - x60*nominal
    # gap, arc span, flat height) and only translates to follow O_p (eccentricity x60) / rotates by
    # phi_k -- so it does not "breathe" and its body sits on O_p, not O_g. The two circles (groove
    # rcut @ O_g, piece r_prime @ O_p) reproduce the true film x60 as their gap.
    half_arc = bush.half_arc_rad()
    curved_gap_mm = (geometry.cutout_radius_m - bush.piece_outer_radius_m) * mmv
    r_prime = rcut - kk * curved_gap_mm
    fig, ax = plt.subplots(figsize=(7.6, 7.2), dpi=140)
    ax.set_aspect("equal")
    ax.grid(True, color="#eef0f3", lw=0.6)
    ogx = kk * groove[0]  # O_g offset from the vane centre-line, exaggerated
    ax.add_patch(plt.Circle((ogx, 0), rcut, fill=False, ec="#c0392b", lw=1.8, zorder=1))
    ax.add_patch(
        plt.Rectangle(
            (-vane_half, -rcut - 1),
            2 * vane_half,
            2 * rcut + 2,
            facecolor=c["vane"],
            ec="#333",
            lw=1,
            zorder=2,
        )
    )
    ax.axvline(0.0, color="#f2c200", ls="--", lw=1.1, zorder=3)  # vane centre-line
    for side, op, cf, ff, phi_k in (
        (1.0, in_op, in_c, in_f, s["in_phi"]),
        (-1.0, out_op, out_c, out_f, s["out_phi"]),
    ):
        cx, cy = kk * op[0], kk * (op[1] - groove[1])
        base = (0.0 if side > 0 else pi) + phi_k
        ts = np.linspace(-half_arc, half_arc, 80)
        arc_x = cx + r_prime * np.cos(base + ts)
        arc_y = cy + r_prime * np.sin(base + ts)
        fx = side * (vane_half + kk * ff * 1e-3)
        ax.fill(
            list(arc_x) + [fx, fx],
            list(arc_y) + [arc_y[-1], arc_y[0]],
            facecolor=c["bush"],
            ec="#222",
            lw=1.2,
            zorder=4,
        )
        ha = "left" if side > 0 else "right"
        ax.annotate(
            f"곡면 {cf:.2f} µm",
            (cx + r_prime * cos(base), cy + r_prime * sin(base)),
            (side * (rcut + 0.5), 4.6),
            fontsize=9,
            color="#c0392b",
            ha=ha,
            arrowprops={"arrowstyle": "->", "color": "#c0392b", "lw": 0.9},
        )
        ax.annotate(
            f"평면 {ff:.2f} µm",
            (fx, cy - 3.0),
            (side * (vane_half + 2.4), -6.8),
            fontsize=9,
            color="#2166ac",
            ha=ha,
            arrowprops={"arrowstyle": "->", "color": "#2166ac", "lw": 0.9},
        )
    ax.plot(ogx, 0, "D", color="#111", ms=7, zorder=6)
    ax.plot(
        kk * in_op[0],
        kk * (in_op[1] - groove[1]),
        "o",
        color="#2166ac",
        ms=7,
        mec="white",
        zorder=6,
    )
    ax.plot(
        kk * out_op[0],
        kk * (out_op[1] - groove[1]),
        "o",
        color="#e8752a",
        ms=7,
        mec="white",
        zorder=6,
    )
    ax.annotate("O_g", (ogx, 0), (ogx - 0.6, 1.6), fontsize=9, ha="right", color="#111")
    ax.annotate(
        "O_p(IN)", (kk * in_op[0], 0), (kk * in_op[0] + 0.5, 2.6), fontsize=9, color="#2166ac"
    )
    ax.annotate(
        "O_p(OUT)",
        (kk * out_op[0], 0),
        (kk * out_op[0] - 0.5, -2.8),
        fontsize=9,
        ha="right",
        color="#e8752a",
    )
    ax.annotate(
        "베인 중선", (0, rcut + 0.3), (0, rcut + 1.3), fontsize=8.5, ha="center", color="#b38f00"
    )
    ax.annotate("베인", (0, -rcut - 0.2), (0, -rcut - 0.9), fontsize=9, color="#eee", ha="center")
    ax.text(
        rcut + 0.6,
        rcut + 1.6,
        f"IN\nO_p-O_g = {s['in_ecc']:.1f} µm\nO_p-베인 = {s['in_vx']:+.0f} µm",
        fontsize=8,
        ha="left",
        va="top",
        color="#333",
        bbox={"boxstyle": "round", "fc": "#fffbe6", "ec": "#d8c98a"},
    )
    ax.text(
        -rcut - 0.6,
        rcut + 1.6,
        f"OUT\nO_p-O_g = {s['out_ecc']:.1f} µm\nO_p-베인 = {s['out_vx']:+.0f} µm",
        fontsize=8,
        ha="right",
        va="top",
        color="#333",
        bbox={"boxstyle": "round", "fc": "#fffbe6", "ec": "#d8c98a"},
    )
    ax.text(
        0,
        -rcut - 2.4,
        f"O_g는 베인 중선에서 {groove[0] * 1e3:+.0f} µm; 중심 치우침 ×{kk} "
        "표시, 몸체 실척·겹침 없음",
        fontsize=8,
        ha="center",
        va="top",
        color="#666",
    )
    ax.set_xlim(-rcut - 4.5, rcut + 4.5)
    ax.set_ylim(-rcut - 3.4, rcut + 3.2)
    ax.set_xlabel("x (mm, 베인 중선 기준)")
    ax.set_title(
        f"θ={s['th']:.0f}° 부시 클리어런스 ×{kk} + 중심 치우침(O_p, O_g)",
        fontsize=12,
        fontweight="bold",
    )
    _save(fig, path)


def render_assembly_journal_clearance(curve: PressureCurve, path: Path) -> None:
    """Crank-pin journal clearance at the most-loaded angle, exaggerated x60 (§4.13).

    The crank pin in the rotor bore, the journal oil film clearance and the pin
    eccentricity exaggerated by the same factor as the bush-clearance figure.
    """

    s = _assembly_state_42(curve)
    Rpin, cj, ej = s["Rpin"], s["cj"], s["ej"]
    c, kk = _ASM_COLORS, _ASM_K
    fig, ax = plt.subplots(figsize=(6.8, 6.8), dpi=140)
    ax.set_aspect("equal")
    ax.grid(True, color="#eef0f3", lw=0.6)
    ax.add_patch(plt.Circle((0, 0), Rpin + kk * cj, fill=False, ec="#999", lw=1.8))
    ax.add_patch(plt.Circle(kk * ej * 1e-3, Rpin, facecolor=c["pin"], ec="#8a7", lw=1.2))
    ax.plot(0, 0, "+", color="#999", ms=9, mew=1.3)
    ax.plot(*(kk * ej * 1e-3), "+", color="#7a5", ms=9, mew=1.3)
    ax.annotate(
        "로터 보어", (0, Rpin + kk * cj), (0, Rpin + kk * cj + 0.8), fontsize=9, ha="center"
    )
    ax.annotate(
        f"크랭크핀 (편심 {np.hypot(*ej):.1f} µm)",
        (0, -Rpin),
        (0, -Rpin - 1.6),
        fontsize=9,
        ha="center",
        color="#7a5",
    )
    ax.annotate(
        f"저널 유막 c_j = {cj * 1e3:.0f} µm",
        (Rpin, 0),
        (Rpin + 0.5, 0),
        fontsize=9,
        va="center",
        color="#333",
    )
    rr = Rpin + kk * cj + 2.2
    ax.set_xlim(-rr, rr)
    ax.set_ylim(-rr, rr)
    ax.set_xlabel("x (mm, 보어 국소)")
    ax.set_title(f"θ={s['th']:.0f}° 저널 클리어런스 ×{kk} 과장", fontsize=12, fontweight="bold")
    _save(fig, path)


def _bush_clearance_history(curve: PressureCurve):
    """Crank angles (sorted, deg) and a per-piece µm-series extractor for the coupled orbit."""

    orbit = _coupled_bush_orbit(curve)
    ang = [degrees(a) for a in orbit.crank_angle_rad]
    order = sorted(range(len(ang)), key=lambda i: ang[i])
    th = [ang[i] for i in order]

    def series(values):
        return [values[i] * 1e6 for i in order]

    return orbit, th, series


def _bush_clearance_axes(ax) -> None:
    ax.grid(True, color="#e2e5ea", lw=0.7)
    ax.set_xlim(0, 360)
    ax.set_xticks(range(0, 361, 90))
    ax.set_xlabel("크랭크각 θ (deg)")
    ax.set_ylabel("최소 유막 두께 h (µm)")


def render_bush_curved_clearance_vs_crank(curve: PressureCurve, path: Path) -> None:
    """Bush curved (rotor-groove) min film thickness of each piece over the cycle (§4.14).

    The min-film *history* the single-angle figure only samples at one instant; the IN and
    OUT pieces alternate as the gas moment reverses, the curved IN film reaching ~0.7 µm.
    """

    geometry = curve.geometry
    bush = SwingBush()
    curved_gap = (geometry.cutout_radius_m - bush.piece_outer_radius_m) * 1e6  # 30 um
    blue, orange, grey = "#2166ac", "#e8752a", "#8a8f98"
    orbit, th, series = _bush_clearance_history(curve)

    fig, ax = plt.subplots(figsize=(8.6, 5.4), dpi=140)
    _bush_clearance_axes(ax)
    for tag, vals, color in (
        ("IN", orbit.in_curved_film_m, blue),
        ("OUT", orbit.out_curved_film_m, orange),
    ):
        y = series(vals)
        ax.plot(th, y, color=color, lw=2.0, label=f"{tag} (min {min(y):.2f} µm)")
        j = min(range(len(y)), key=lambda i: y[i])
        ax.plot(th[j], y[j], "v", color=color, ms=8)
    ax.axhline(curved_gap, ls="--", color=grey, lw=1.3, label=f"공칭 gap {curved_gap:.0f} µm")
    ax.set_title(
        "부시 곡면막 (로터 그루브) 클리어런스 이력 — 사이클 전체", fontsize=12, fontweight="bold"
    )
    ax.legend(fontsize=9)
    _save(fig, path)


def render_bush_flat_clearance_vs_crank(curve: PressureCurve, path: Path) -> None:
    """Bush flat (vane) min film thickness of each piece over the cycle (§4.14)."""

    geometry = curve.geometry
    bush = SwingBush()
    flat_gap = film_thicknesses_m(geometry, bush)[0] * 1e6  # 10 um
    blue, orange, grey = "#2166ac", "#e8752a", "#8a8f98"
    orbit, th, series = _bush_clearance_history(curve)

    fig, ax = plt.subplots(figsize=(8.6, 5.4), dpi=140)
    _bush_clearance_axes(ax)
    for tag, vals, color in (
        ("IN", orbit.in_flat_film_m, blue),
        ("OUT", orbit.out_flat_film_m, orange),
    ):
        y = series(vals)
        ax.plot(th, y, color=color, lw=2.0, label=f"{tag} (min {min(y):.2f} µm)")
        j = min(range(len(y)), key=lambda i: y[i])
        ax.plot(th[j], y[j], "v", color=color, ms=8)
    ax.axhline(flat_gap, ls="--", color=grey, lw=1.3, label=f"공칭 c_f {flat_gap:.0f} µm")
    ax.set_title(
        "부시 평면막 (베인) 클리어런스 이력 — 사이클 전체", fontsize=12, fontweight="bold"
    )
    ax.legend(fontsize=9)
    _save(fig, path)


# ==========================================================================
# Raw-data (.dat) export for research post-processing (Tecplot ASCII).
#
# The PNGs above are illustration; these functions dump the numbers behind the
# coupled-orbit and validation figures. Two shapes: 1-D point zones (a quantity
# vs crank angle / eccentricity) and 2-D ordered surface zones (film thickness
# over crank-angle x position). Registered in DATASETS, emitted with --data.
# ==========================================================================


def _orbit_common(curve: PressureCurve) -> dict:
    """Coupled-orbit handle plus the fixed bush/journal geometry the exporters need."""

    geometry = curve.geometry
    orbit = _coupled_bush_orbit(curve)
    bush = SwingBush()
    return {
        "geometry": geometry,
        "orbit": orbit,
        "bush": bush,
        "half_arc": bush.half_arc_rad(),
        "curved_gap": geometry.cutout_radius_m - bush.piece_outer_radius_m,
        "flat_gap": film_thicknesses_m(geometry, bush)[0],
        "shift": bush.piece_shift_m,
        "lever": geometry.cutout_offset_m,
        "throw": geometry.eccentricity_m,
        "c_j": orbit.journal_clearance_m,
    }


def _orbit_frames(c: dict):
    """Per-sample crank angle, rotor eccentricity, and groove centre (world), from the orbit.

    Returns arrays sorted by crank angle so the 2-D field zones have a monotonic J-axis.
    """

    geometry, orbit = c["geometry"], c["orbit"]
    theta = np.array(orbit.crank_angle_rad)
    order = np.argsort(theta)
    theta = theta[order]
    ejx = np.array(orbit.eccentricity_x_m)[order]
    ejy = np.array(orbit.eccentricity_y_m)[order]
    phi = (
        np.array([prescribed_state(geometry, float(t)).rotor_orientation_rad for t in theta])
        + np.array(orbit.rotor_attitude_deviation_rad)[order]
    )
    bore_x = c["throw"] * np.sin(theta) + ejx
    bore_y = c["throw"] * np.cos(theta) + ejy
    groove_x = bore_x + c["lever"] * np.cos(phi)
    groove_y = bore_y + c["lever"] * np.sin(phi)
    return {"order": order, "theta": theta, "ejx": ejx, "ejy": ejy, "gx": groove_x, "gy": groove_y}


def _write_film_strand(path, x_name, frames, title):
    """Write one 1-D line zone per crank angle as a Tecplot transient *strand*.

    Each frame is ``(solution_time_s, angle_deg, [x_column, h_um_column])``; every zone shares
    one STRANDID and carries its own SOLUTIONTIME, so Tecplot animates through crank angle
    (x = position, y = h_um). Frames with non-increasing solution time (the wrap seam) are
    skipped so the strand is strictly ordered.
    """

    # SOLUTIONTIME is the CRANK ANGLE in degrees (not seconds): Tecplot's time slider then
    # reads the crank angle directly (0..360). Drop frames that collide on the written value
    # (the revolution wrap seam samples the same crank angle twice, distinct only at ~1e-14),
    # keying on the printed representation so the strand is strictly increasing.
    zones = []
    seen: set[str] = set()
    for _time_s, angle_deg, columns in frames:
        key = f"{angle_deg:.9g}"
        if key in seen:
            continue
        seen.add(key)
        zones.append(
            point_zone(f"theta_{angle_deg:.1f}deg", columns, strand_id=1, solution_time=angle_deg)
        )
    write_dat(path, [x_name, "h_um"], zones, title=title)


def data_orbit_film_journal(curve: PressureCurve, path: Path) -> None:
    """Transient: journal film h(phi) as one line-zone per crank angle (step/animate in Tecplot)."""

    c = _orbit_common(curve)
    fr = _orbit_frames(c)
    time_s = fr["theta"] / c["geometry"].angular_speed_rad_s
    angle = np.degrees(fr["theta"])
    phi_deg = np.degrees(np.linspace(0.0, 2.0 * pi, 361))
    phi = np.radians(phi_deg)
    ex, ey = fr["ejx"][None, :], fr["ejy"][None, :]
    h_um = (c["c_j"] - (ex * np.cos(phi)[:, None] + ey * np.sin(phi)[:, None])) * 1e6  # [phi, time]
    frames = [
        (float(time_s[j]), float(angle[j]), [phi_deg, h_um[:, j]]) for j in range(time_s.size)
    ]
    _write_film_strand(path, "phi_deg", frames, "Journal film h(phi) per crank angle (transient)")


def _data_curved_map(curve, side, x_key, y_key, phi_key, path, title):
    c = _orbit_common(curve)
    fr = _orbit_frames(c)
    order = fr["order"]
    orbit = c["orbit"]
    time_s = fr["theta"] / c["geometry"].angular_speed_rad_s
    angle = np.degrees(fr["theta"])
    xk = np.array(getattr(orbit, x_key))[order]
    yk = np.array(getattr(orbit, y_key))[order]
    phik = np.array(getattr(orbit, phi_key))[order]
    ecc_x, ecc_y = xk - fr["gx"], yk - fr["gy"]
    center = (0.0 if side > 0.0 else pi) + phik  # arc centre azimuth per sample
    beta_rel_deg = np.degrees(np.linspace(-c["half_arc"], c["half_arc"], 181))
    beta = center[None, :] + np.radians(beta_rel_deg)[:, None]  # [beta, time]
    h_um = (c["curved_gap"] - (ecc_x[None, :] * np.cos(beta) + ecc_y[None, :] * np.sin(beta))) * 1e6
    frames = [
        (float(time_s[j]), float(angle[j]), [beta_rel_deg, h_um[:, j]]) for j in range(time_s.size)
    ]
    _write_film_strand(path, "beta_rel_deg", frames, title)


def data_orbit_film_curved_in(curve: PressureCurve, path: Path) -> None:
    """Transient: IN-piece curved (rotor-groove) film h(beta) per crank angle."""
    _data_curved_map(
        curve, 1.0, "in_piece_x_m", "in_piece_y_m", "in_piece_attitude_rad", path,
        "IN bush curved film h(beta) per crank angle (transient)",
    )


def data_orbit_film_curved_out(curve: PressureCurve, path: Path) -> None:
    """Transient: OUT-piece curved (rotor-groove) film h(beta) per crank angle."""
    _data_curved_map(
        curve, -1.0, "out_piece_x_m", "out_piece_y_m", "out_piece_attitude_rad", path,
        "OUT bush curved film h(beta) per crank angle (transient)",
    )


def _data_flat_map(curve, side, x_key, phi_key, path, title):
    c = _orbit_common(curve)
    fr = _orbit_frames(c)
    order = fr["order"]
    geometry, orbit, bush = c["geometry"], c["orbit"], c["bush"]
    time_s = fr["theta"] / geometry.angular_speed_rad_s
    angle = np.degrees(fr["theta"])
    xk = np.array(getattr(orbit, x_key))[order]
    phik = np.array(getattr(orbit, phi_key))[order]
    lf = np.array([flat_contact_length_m(geometry, float(t), bush) for t in fr["theta"]])
    xi = np.linspace(-0.5, 0.5, 101)  # pad fraction; the actual mm position (s) varies with θ
    approach = c["shift"] - side * xk
    frames = []
    for j in range(time_s.size):
        s_mm = xi * lf[j] * 1e3  # actual pad position at this crank angle
        h_um = ((c["flat_gap"] - approach[j]) + xi * lf[j] * phik[j]) * 1e6
        frames.append((float(time_s[j]), float(angle[j]), [s_mm, h_um]))
    _write_film_strand(path, "s_mm", frames, title)


def data_orbit_film_flat_in(curve: PressureCurve, path: Path) -> None:
    """Transient: IN-piece flat (vane) film h(s) per crank angle."""
    _data_flat_map(
        curve, 1.0, "in_piece_x_m", "in_piece_attitude_rad", path,
        "IN bush flat film h(s) per crank angle (transient)",
    )


def data_orbit_film_flat_out(curve: PressureCurve, path: Path) -> None:
    """Transient: OUT-piece flat (vane) film h(s) per crank angle."""
    _data_flat_map(
        curve, -1.0, "out_piece_x_m", "out_piece_attitude_rad", path,
        "OUT bush flat film h(s) per crank angle (transient)",
    )


def data_orbit_state_timeseries(curve: PressureCurve, path: Path) -> None:
    """1-D series: rotor/piece kinematics, film thicknesses, and seal contact vs crank angle."""

    orbit = _coupled_bush_orbit(curve)
    theta = np.degrees(np.array(orbit.crank_angle_rad))
    order = np.argsort(theta)
    n = theta.size

    def col(seq, scale=1.0):
        arr = np.array(seq, dtype=float)
        return (arr[order] * scale) if arr.size == n else np.zeros(n)

    columns = [
        theta[order],
        col(orbit.eccentricity_x_m, 1e6), col(orbit.eccentricity_y_m, 1e6),
        col(orbit.rotor_attitude_deviation_rad, 1e3),
        col(orbit.in_piece_x_m, 1e3), col(orbit.in_piece_y_m, 1e3),
        col(orbit.in_piece_attitude_rad, 1e3),
        col(orbit.out_piece_x_m, 1e3), col(orbit.out_piece_y_m, 1e3),
        col(orbit.out_piece_attitude_rad, 1e3),
        col(orbit.in_curved_film_m, 1e6), col(orbit.in_flat_film_m, 1e6),
        col(orbit.out_curved_film_m, 1e6), col(orbit.out_flat_film_m, 1e6),
        col(orbit.seal_normal_force_n), col(orbit.seal_penetration_m, 1e6),
    ]
    variables = [
        "theta_deg", "e_jx_um", "e_jy_um", "dphi_r_mrad",
        "in_x_mm", "in_y_mm", "in_phi_mrad", "out_x_mm", "out_y_mm", "out_phi_mrad",
        "in_curved_um", "in_flat_um", "out_curved_um", "out_flat_um",
        "seal_normal_n", "seal_penetration_um",
    ]
    write_dat(path, variables, [point_zone("orbit_state", columns)],
              title="Coupled rotor-bush orbit - kinematics, films, seal (final revolution)")


def data_orbit_force_timeseries(curve: PressureCurve, path: Path) -> None:
    """1-D series: every per-part force/moment channel vs crank angle (N, N*m)."""

    orbit = _coupled_bush_orbit(curve)
    channels = orbit.sample_channels or {}
    theta = np.degrees(np.array(orbit.crank_angle_rad))
    order = np.argsort(theta)
    names = list(channels)
    columns = [theta[order]] + [np.array(channels[k], dtype=float)[order] for k in names]
    variables = ["theta_deg", *names]
    write_dat(path, variables, [point_zone("orbit_forces", columns)],
              title="Coupled rotor-bush orbit - per-part forces and moments (final revolution)")


def data_reynolds_curved(geometry: RotaryGeometry, path: Path) -> None:
    """1-D series: bush curved film — arc_film vs long-bearing load over ε."""
    d = _reynolds_curved_data(geometry)
    columns = [d["eps"], d["lb_mag"], d["arc_mag"]]
    write_dat(path, ["eps", "F_long_bearing_N", "F_arc_film_N"], [point_zone("curved", columns)],
              title=f"Bush curved film validation (max rel. err {d['err']:.2e})")


def data_reynolds_flat(geometry: RotaryGeometry, path: Path) -> None:
    """1-D series: bush flat film — slider_film vs fixed-incline slider over the wedge tilt."""
    d = _reynolds_flat_data(geometry)
    columns = [[t * 1e3 for t in d["tilts"]], d["inc"], d["sld"]]
    write_dat(path, ["tilt_mrad", "W_incline_N", "W_slider_N"], [point_zone("flat", columns)],
              title=f"Bush flat film validation (max rel. err {d['err']:.2e})")


def data_reynolds_journal(geometry: RotaryGeometry, path: Path) -> None:
    """1-D series: journal film — 1-D numerical Reynolds vs Ocvirk short bearing over ε."""
    d = _reynolds_journal_data()
    columns = [d["epsj"], d["ocv"], d["num"]]
    write_dat(path, ["eps", "F_ocvirk_N", "F_reynolds_1d_N"], [point_zone("journal", columns)],
              title=f"Journal film validation (max rel. err {d['err']:.2e})")


# Raw-data manifest (paths are relative to results/data/). Same kind dispatch as FIGURES:
# "curve" reuses the shared coupled orbit, "geometry" needs only the static geometry.
DATASETS: dict[str, tuple[str, object]] = {
    "orbit/state_timeseries.dat": ("curve", data_orbit_state_timeseries),
    "orbit/force_timeseries.dat": ("curve", data_orbit_force_timeseries),
    "orbit/film_journal.dat": ("curve", data_orbit_film_journal),
    "orbit/film_curved_in.dat": ("curve", data_orbit_film_curved_in),
    "orbit/film_curved_out.dat": ("curve", data_orbit_film_curved_out),
    "orbit/film_flat_in.dat": ("curve", data_orbit_film_flat_in),
    "orbit/film_flat_out.dat": ("curve", data_orbit_film_flat_out),
    "validation/reynolds_curved.dat": ("geometry", data_reynolds_curved),
    "validation/reynolds_flat.dat": ("geometry", data_reynolds_flat),
    "validation/reynolds_journal.dat": ("geometry", data_reynolds_journal),
}


FIGURES: dict[str, tuple[str, object]] = {
    # Physics figures — one single-axes figure per file, grouped by topic folder.
    "chamber_pressure/vs_crank.png": ("curve", render_chamber_pressures),
    "indicated_work/pv_diagram.png": ("curve", render_indicated_pv_diagram),
    "indicated_work/thermo_crosscheck.png": ("curve", render_thermo_crosscheck),
    "bush_film/film_pressure.png": ("curve", render_bush_film_pressure),
    "bush_film/sliding_velocity.png": ("curve", render_bush_film_velocity),
    "bush_film/shear_traction.png": ("curve", render_bush_film_shear),
    "bush_film/uniform_clearance_model.png": ("geometry", render_bush_uniform_clearance_model),
    "gas_force/force_components.png": ("curve", render_gas_force_components),
    "gas_force/force_magnitude.png": ("curve", render_gas_force_magnitude),
    "gas_force/force_magnitude_true.png": ("curve", render_gas_force_magnitude_true),
    "gas_force/force_magnitude_circular.png": ("curve", render_gas_force_magnitude_circular),
    "gas_force/gas_torque.png": ("curve", render_gas_torque),
    "gas_force/gas_torque_true.png": ("curve", render_gas_torque_true),
    "gas_force/gas_torque_circular.png": ("curve", render_gas_torque_circular),
    "gas_force/free_body.png": ("curve", render_force_diagram),
    "bearing_load/journal_load_polar.png": ("curve", render_journal_load_polar),
    "bearing_load/journal_load_polar_true.png": ("curve", render_journal_load_polar_true),
    "bearing_load/journal_load_polar_circular.png": ("curve", render_journal_load_polar_circular),
    "bearing_load/drive_torque.png": ("curve", render_drive_torque),
    "bearing_load/drive_torque_circular.png": ("curve", render_drive_torque_circular),
    "bearing_load/friction_breakdown_true.png": ("curve", render_friction_breakdown_true),
    "bearing_load/friction_breakdown_circular.png": ("curve", render_friction_breakdown_circular),
    "bearing_load/journal_concentric_clearance_model.png": (
        "geometry",
        render_journal_concentric_clearance_model,
    ),
    "bearing_load/journal_relative_speed.png": ("geometry", render_journal_relative_speed),
    "bearing_load/journal_eccentricity.png": ("curve", render_journal_eccentricity),
    "bearing_load/reynolds_1d_validation.png": ("geometry", render_reynolds_1d_validation),
    "bearing_load/rotor_orbit.png": ("curve", render_rotor_orbit),
    "bearing_load/eccentric_friction_power.png": ("curve", render_eccentric_friction_power),
    "bush_film/rotor_swing_vs_bush.png": ("geometry", render_rotor_swing_vs_bush),
    "bush_film/bush_translation_vs_vane.png": ("geometry", render_bush_translation_vs_vane),
    "bush_film/reynolds_curved_vs_long_bearing.png": ("geometry", render_film_reynolds_curved),
    "bush_film/reynolds_flat_vs_incline_slider.png": ("geometry", render_film_reynolds_flat),
    "bush_film/reynolds_journal_vs_ocvirk.png": ("geometry", render_film_reynolds_journal),
    "bush_film/film_clearance_journal.png": ("curve", render_film_clearance_journal),
    "bush_film/film_clearance_bush_curved.png": ("curve", render_film_clearance_bush_curved),
    "bush_film/film_clearance_bush_flat.png": ("curve", render_film_clearance_bush_flat),
    "bush_film/bush_curved_clearance_vs_crank.png": (
        "curve",
        render_bush_curved_clearance_vs_crank,
    ),
    "bush_film/bush_flat_clearance_vs_crank.png": ("curve", render_bush_flat_clearance_vs_crank),
    "assembly/layout_42deg.png": ("curve", render_assembly_layout),
    "assembly/bush_clearance_42deg.png": ("curve", render_assembly_bush_clearance),
    "assembly/journal_clearance_42deg.png": ("curve", render_assembly_journal_clearance),
    "bearing_load/friction_dynamic_vs_quasistatic.png": (
        "curve",
        render_friction_dynamic_vs_quasistatic,
    ),
    "bearing_load/roller_vs_crank.png": ("geometry", render_roller_vs_crank),
    "volumetric_efficiency/recompression_pressure.png": ("curve", render_recompression_pressure),
    "volumetric_efficiency/chamber_mass.png": ("curve", render_chamber_mass),
    "reed_valve/discharge_pressure.png": ("curve", render_reed_valve_pressure),
    "reed_valve/port_open_area.png": ("curve", render_reed_valve_port_area),
    "rotor_motion/mechanism.gif": ("curve", render_rotor_motion),
    # Geometry drawings — dimensioned montages kept whole, under one folder.
    "geometry/port_geometry.png": ("geometry", render_port_geometry),
    "geometry/tdc_bdc_definition.png": ("geometry", render_tdc_bdc_definition),
    "geometry/crank_angle_definition.png": ("geometry", render_crank_angle_definition),
    "geometry/journal_film_coordinates.png": ("geometry", render_journal_film_coordinates),
    "geometry/journal_film_axial.png": ("geometry", render_journal_film_axial),
    "geometry/bush_film_coordinates.png": ("geometry", render_bush_film_coordinates),
    "geometry/bush_attitude_reference.png": ("geometry", render_bush_attitude_reference),
    "geometry/bush_2piece_coordinates.png": ("geometry", render_bush_2piece_coordinates),
    "geometry/rotor_mouth_lip_detail.png": ("geometry", render_rotor_mouth_lip_detail),
    "geometry/dimensioned_side_section.png": ("geometry", render_dimensioned_side_section),
    "geometry/stepped_vane_structure.png": ("geometry", render_stepped_vane_structure),
    "geometry/bush_placement_clearances.png": ("geometry", render_bush_placement_clearances),
    "geometry/dimensioned_top_view.png": ("geometry", render_dimensioned_top_view),
    "geometry/vane_side_view.png": ("geometry", render_vane_side_view),
    "geometry/vane_model_comparison.png": ("geometry", render_vane_model_comparison),
    "geometry/chamber_case1_seal_over.png": ("geometry", render_chamber_case1_seal_over),
    "geometry/chamber_case2_gap.png": ("geometry", render_chamber_case2_gap),
    "geometry/geometry_master.png": ("geometry", render_geometry_master),
}

# Legacy figures still awaiting a renderer, kept until ported. Now empty: every
# legacy figure has been ported into FIGURES, except rotor_mouth_48deg.png, which
# drew the superseded SYMMETRIC 48°/48° mouth. Its premise is obsolete (the mouth
# is asymmetric), so it is deliberately left outside the manifest for --prune to
# reclaim; the current mouth is shown by geometry/port_geometry.png and
# geometry/rotor_mouth_lip_detail.png.
PENDING_LEGACY_FIGURES: tuple[str, ...] = ()

# The set --prune keeps: figures this script produces plus the not-yet-ported
# legacy figures it deliberately preserves.
GENERATED_FIGURES: tuple[str, ...] = tuple(FIGURES) + PENDING_LEGACY_FIGURES


def prune_stale(prune: bool) -> list[Path]:
    """List, and optionally delete, images (any subfolder) this script does not produce."""

    kept = {name.lower() for name in GENERATED_FIGURES}
    suffixes = {".png", ".gif", ".jpg", ".jpeg", ".svg", ".mp4"}
    stale = [
        item
        for item in sorted(RESULTS_DIR.rglob("*"))
        if item.is_file()
        and item.suffix.lower() in suffixes
        and item.relative_to(RESULTS_DIR).as_posix().lower() not in kept
    ]
    for item in stale:
        rel = item.relative_to(RESULTS_DIR).as_posix()
        if prune:
            item.unlink()
            print(f"  deleted stale figure: {rel}")
        else:
            print(f"  stale (run with --prune to delete): {rel}")
    if prune:
        # Remove any directories left empty after pruning (deepest first).
        for item in sorted(RESULTS_DIR.rglob("*"), reverse=True):
            if item.is_dir():
                try:
                    item.rmdir()
                except OSError:
                    pass
    return stale


def _ensure_curve(curve: PressureCurve | None, geometry: RotaryGeometry) -> PressureCurve:
    """Build the shared port-timed pressure trace once, reusing it across figures and data."""

    if curve is None:
        print("Building the port-timed cycle trace (a minute or so)...")
        curve = PressureCurve(geometry)
    return curve


def _run_manifest(manifest: dict, base_dir: Path, geometry, curve, verb: str):
    """Run a FIGURES/DATASETS manifest, dispatching by kind; returns the (reused) curve."""

    for name, (kind, fn) in manifest.items():
        path = base_dir / name
        path.parent.mkdir(parents=True, exist_ok=True)
        if kind == "curve":
            curve = _ensure_curve(curve, geometry)
            print(f"{verb} {name} ...")
            fn(curve, path)
        else:
            print(f"{verb} {name} ...")
            fn(geometry, path)
    return curve


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--prune",
        action="store_true",
        help="delete results/ images this script does not produce",
    )
    parser.add_argument(
        "--data",
        action="store_true",
        help="also export raw Tecplot .dat data under results/data/ (illustration PNGs + data)",
    )
    parser.add_argument(
        "--data-only",
        action="store_true",
        help="export only the raw .dat data (skip figures and pruning)",
    )
    args = parser.parse_args()

    RESULTS_DIR.mkdir(exist_ok=True)
    _use_korean_font()
    geometry = RotaryGeometry.default()

    curve: PressureCurve | None = None
    if not args.data_only:
        curve = _run_manifest(FIGURES, RESULTS_DIR, geometry, curve, "Rendering")

    if args.data or args.data_only:
        # Reuses the same curve (and its cached coupled orbit) as the figures above.
        curve = _run_manifest(DATASETS, DATA_DIR, geometry, curve, "Writing data")

    if not args.data_only:
        if PENDING_LEGACY_FIGURES:
            print(f"Preserving {len(PENDING_LEGACY_FIGURES)} legacy figure(s) not yet ported.")
        print("Checking for figures outside the manifest...")
        prune_stale(args.prune)

    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
