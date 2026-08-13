"""Bush-focused animation: the rotor mouth (real rotor shape) with a x60 clearance inset.

A **separate** viewer from :mod:`mochi.gui` (the prescribed-motion test GUI). This one
animates the **coupled 9-DOF rotor+bush orbit** over one revolution and shows, for every
crank angle:

* **Main panel (zoomed out, true scale):** the real rotor contour at the mouth
  (:func:`mochi.rotor_profile.rotor_contour`, whose groove *is* part of the rotor shape),
  the fixed vane, and the two swing-bush pieces -- so the rotor behaviour near the bush is
  visible as the machine turns. Clearances are micrometre-scale and invisible here by design.
* **Inset (x60 exaggerated):** the same tight, vane-centred view as the
  ``assembly/bush_clearance`` figure -- the two pieces' curved (rotor-groove) and flat (vane)
  films and the centre offsets ``O_g``, ``O_p(IN/OUT)`` blown up by a uniform factor so the
  clearance behaviour is legible.

The orbit is expensive (~minutes); it is computed once (optionally cached to ``.npz``) and the
GUI/animation then scrubs the cached frames. Two entry points share one drawing routine:

* :func:`render_bush_cycle_gif` -- write an animated GIF (matplotlib ``PillowWriter``).
* :class:`BushViewerApp` / :func:`main` -- a Tk window embedding the same figure, with a
  crank-angle slider, play/pause, and a *Save GIF* button.
"""

from __future__ import annotations

import argparse
from collections.abc import Callable
from dataclasses import dataclass
from math import cos, pi, sin
from pathlib import Path
from typing import TYPE_CHECKING, Any

import matplotlib
import numpy as np

from mochi.bush_outline import bush_piece_outline_mm, vane_outline_mm
from mochi.chamber_volume import SwingBush
from mochi.chambers import build_cycle_trace
from mochi.kinematics import MM, RotaryGeometry, prescribed_state
from mochi.rotor_bush_dynamics import integrate_rotor_bush_orbit
from mochi.rotor_profile import rotor_contour

if TYPE_CHECKING:
    import tkinter as tk

    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

# Uniform clearance/centre-offset exaggeration in the inset, matching the static
# ``assembly/bush_clearance`` figure (``generate_results._ASM_K``). Held fixed for every
# crank angle so frames are comparable.
EXAGGERATION = 60

_COLORS = {
    "rotor": "#d9d9d9",
    "rotor_edge": "#555555",
    "bush": "#aeb8c7",
    "vane": "#5a5a5a",
    "groove": "#c0392b",
    "in": "#2166ac",
    "out": "#e8752a",
    "guide": "#8a8f98",
}


def _set_korean_font() -> None:
    """Prefer a Hangul-capable font so the Korean labels render (falls back to DejaVu)."""

    import matplotlib.pyplot as plt
    from matplotlib import font_manager

    installed = {f.name for f in font_manager.fontManager.ttflist}
    for name in ("Malgun Gothic", "Gulim", "NanumGothic", "AppleGothic"):
        if name in installed:
            plt.rcParams["font.family"] = [name, "DejaVu Sans"]
            break
    plt.rcParams["axes.unicode_minus"] = False


# ---------------------------------------------------------------------------
# Orbit data (compute once, optionally cached)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class BushOrbitData:
    """The per-frame coupled-orbit fields the animation needs (SI units)."""

    crank: np.ndarray
    ejx: np.ndarray
    ejy: np.ndarray
    rdev: np.ndarray
    inx: np.ndarray
    iny: np.ndarray
    inphi: np.ndarray
    outx: np.ndarray
    outy: np.ndarray
    outphi: np.ndarray
    inc: np.ndarray
    inf: np.ndarray
    outc: np.ndarray
    outf: np.ndarray
    cj: float

    @property
    def frames(self) -> int:
        return int(self.crank.size)


def compute_orbit_data(
    geometry: RotaryGeometry,
    *,
    revolutions: int = 4,
    samples: int = 180,
    cache_path: str | Path | None = None,
    recompute: bool = False,
) -> BushOrbitData:
    """Coupled 9-DOF orbit fields, loaded from ``cache_path`` if present else computed.

    The coupled integration takes minutes; pass ``cache_path`` to persist the extracted
    arrays as an ``.npz`` and reuse them on later launches (``recompute=True`` forces a
    fresh integration and overwrites the cache).
    """

    if cache_path is not None:
        cache_path = Path(cache_path)
        if cache_path.exists() and not recompute:
            data = np.load(cache_path)
            return BushOrbitData(
                **{k: data[k] for k in data.files if k != "cj"}, cj=float(data["cj"])
            )

    trace = build_cycle_trace(geometry)
    orbit = integrate_rotor_bush_orbit(
        geometry,
        revolutions=revolutions,
        samples=samples,
        grid_samples=180,
        n_beta=121,
        n_s=81,
        trace=trace,
        seal_contact=True,
    )
    result = BushOrbitData(
        crank=np.array(orbit.crank_angle_rad),
        ejx=np.array(orbit.eccentricity_x_m),
        ejy=np.array(orbit.eccentricity_y_m),
        rdev=np.array(orbit.rotor_attitude_deviation_rad),
        inx=np.array(orbit.in_piece_x_m),
        iny=np.array(orbit.in_piece_y_m),
        inphi=np.array(orbit.in_piece_attitude_rad),
        outx=np.array(orbit.out_piece_x_m),
        outy=np.array(orbit.out_piece_y_m),
        outphi=np.array(orbit.out_piece_attitude_rad),
        inc=np.array(orbit.in_curved_film_m),
        inf=np.array(orbit.in_flat_film_m),
        outc=np.array(orbit.out_curved_film_m),
        outf=np.array(orbit.out_flat_film_m),
        cj=float(orbit.journal_clearance_m),
    )
    if cache_path is not None:
        np.savez(
            cache_path,
            **{f.name: getattr(result, f.name) for f in result.__dataclass_fields__.values()},
        )
    return result


# ---------------------------------------------------------------------------
# Per-frame geometric state
#
# The bush-piece and vane outlines (mm) are shared with the static ``results/`` figures via
# :mod:`mochi.bush_outline`. This viewer draws the vane with ``close_bore_arc=False`` so the
# zoomed-out tongue closes with a straight chord, as it always has.
# ---------------------------------------------------------------------------


def _frame_state(geometry: RotaryGeometry, data: BushOrbitData, j: int) -> dict[str, Any]:
    """Lab-frame positions (mm) and clearances (um) for frame ``j`` of the orbit."""

    mm = 1.0 / MM
    theta = float(data.crank[j])
    state = prescribed_state(geometry, theta)
    phi = state.rotor_orientation_rad + float(data.rdev[j])
    throw = geometry.eccentricity_m
    pin = np.array([throw * sin(theta), throw * cos(theta)])
    bore = pin + np.array([data.ejx[j], data.ejy[j]])
    groove = bore + geometry.cutout_offset_m * np.array([cos(phi), sin(phi)])
    return {
        "theta_deg": (theta * 180.0 / pi) % 360.0,
        "state": state,
        "groove": groove * mm,
        "in_op": np.array([data.inx[j], data.iny[j]]) * mm,
        "out_op": np.array([data.outx[j], data.outy[j]]) * mm,
        "in_c": float(data.inc[j]) * 1e6,
        "in_f": float(data.inf[j]) * 1e6,
        "out_c": float(data.outc[j]) * 1e6,
        "out_f": float(data.outf[j]) * 1e6,
        "in_phi": float(data.inphi[j]),
        "out_phi": float(data.outphi[j]),
        "vane_tip_y_m": state.vane_tip_m[1],
    }


# ---------------------------------------------------------------------------
# Drawing (one crank angle onto the two axes)
# ---------------------------------------------------------------------------


def _main_window_mm(
    geometry: RotaryGeometry, data: BushOrbitData
) -> tuple[float, float, float, float]:
    """Fixed zoomed-out window (mm) framing the rotor mouth across the whole cycle."""

    mm = 1.0 / MM
    gy = []
    for j in range(data.frames):
        s = _frame_state(geometry, data, j)
        gy.append(s["groove"][1])
    rcut = geometry.cutout_radius_m * mm
    y_lo = min(gy) - rcut - 5.0
    y_hi = geometry.cylinder_radius_m * mm + 3.0
    half_w = 0.5 * (y_hi - y_lo) * 1.02  # keep the aspect near square
    return -half_w, half_w, y_lo, y_hi


def draw_main(
    ax: Axes,
    geometry: RotaryGeometry,
    data: BushOrbitData,
    j: int,
    window: tuple[float, float, float, float],
) -> None:
    """Zoomed-out, true-scale lab view: real rotor shape at the mouth, vane, bush pieces."""

    from matplotlib.patches import Circle

    mm = 1.0 / MM
    s = _frame_state(geometry, data, j)
    theta_rad = float(data.crank[j])
    c = _COLORS
    ax.clear()
    ax.set_aspect("equal")
    ax.grid(True, color="#eef0f3", lw=0.6)

    # Cylinder wall (context; only the top arc lands in the window).
    ax.add_patch(Circle((0, 0), geometry.cylinder_radius_m * mm, fill=False, ec="#333", lw=1.4))

    # Real rotor contour -- the groove pocket is part of this shape.
    contour = rotor_contour(geometry, theta_rad)
    mat = np.array(contour.material) * mm
    ax.fill(mat[:, 0], mat[:, 1], facecolor=c["rotor"], ec="none", zorder=1)
    for edge in (contour.od_arc, contour.inlet_flat, contour.mouth_path):
        e = np.array(edge) * mm
        ax.plot(e[:, 0], e[:, 1], color=c["rotor_edge"], lw=1.6, zorder=2, solid_joinstyle="round")

    # Fixed vane.
    vx, vy = vane_outline_mm(geometry, s["vane_tip_y_m"], close_bore_arc=False)
    ax.fill(vx, vy, facecolor=c["vane"], ec="#333", lw=1.0, zorder=4)

    # Two swing-bush pieces at the prescribed cutout centre (true scale, aligned with the
    # rotor contour above; the micrometre orbit offset is the inset's job).
    groove_x_m, groove_y_m = s["state"].cutout_center_m
    bush = SwingBush()
    for side in (1.0, -1.0):
        bx, by = bush_piece_outline_mm(groove_x_m + side * bush.piece_shift_m, groove_y_m, side)
        ax.fill(bx, by, facecolor=c["bush"], ec="#222", lw=1.0, zorder=5)

    # Centre marker O_g (groove centre).
    ox_mm, oy_mm = groove_x_m / MM, groove_y_m / MM
    ax.plot(ox_mm, oy_mm, "D", color=c["groove"], ms=6, mec="white", zorder=6)
    ax.annotate("O_g", (ox_mm, oy_mm), (ox_mm + 0.8, oy_mm + 1.2), fontsize=8.5, color=c["groove"])

    x0, x1, y0, y1 = window
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_xlabel("x (mm, 실린더 좌표)")
    ax.set_ylabel("y (mm)")
    ax.set_title(
        f"로터 마우스 실척 — θ = {s['theta_deg']:.0f}°  (로터·부시 실제 거동)",
        fontsize=11,
        fontweight="bold",
    )


def draw_inset(ax: Axes, geometry: RotaryGeometry, data: BushOrbitData, j: int) -> None:
    """Tight, vane-centred, x60 view: curved/flat films and O_g/O_p offsets (bush_clearance)."""

    from matplotlib.patches import Circle, Rectangle

    s = _frame_state(geometry, data, j)
    bush = SwingBush()
    mm = 1.0 / MM
    rcut = geometry.cutout_radius_m * mm
    vane_half = 0.5 * geometry.vane_width_m * mm
    groove = s["groove"]
    c, kk = _COLORS, EXAGGERATION

    # Rigid caricature: the piece keeps a CONSTANT shape (outer radius, arc span, flat height)
    # and only translates to follow O_p (eccentricity x60) / rotates by phi_k. Shrinking the
    # drawn outer radius by the exaggerated *nominal* gap makes the film a visible band while the
    # two circles (groove rcut @ O_g, piece r_prime @ O_p) reproduce the true film x60 as their
    # gap. Contrast the old schematic, which set the radius from the live film -> the piece
    # appeared to breathe and sat on O_g instead of O_p.
    curved_gap_mm = (geometry.cutout_radius_m - bush.piece_outer_radius_m) * mm
    half_arc = bush.half_arc_rad()
    r_prime = rcut - kk * curved_gap_mm  # constant caricature outer radius
    ax.clear()
    ax.set_aspect("equal")
    ax.grid(True, color="#eef0f3", lw=0.6)

    ogx = kk * groove[0]  # O_g offset from the vane centre-line, exaggerated
    ax.add_patch(Circle((ogx, 0), rcut, fill=False, ec=c["groove"], lw=1.6, zorder=1))
    ax.add_patch(
        Rectangle(
            (-vane_half, -rcut - 1),
            2 * vane_half,
            2 * rcut + 2,
            facecolor=c["vane"],
            ec="#333",
            lw=1,
            zorder=2,
        )
    )
    ax.axvline(0.0, color="#f2c200", ls="--", lw=1.0, zorder=3)  # vane centre-line
    pieces = (
        (1.0, s["in_op"], s["in_c"], s["in_f"], s["in_phi"], c["in"], "O_p(IN)"),
        (-1.0, s["out_op"], s["out_c"], s["out_f"], s["out_phi"], c["out"], "O_p(OUT)"),
    )
    for side, op, cf, ff, phi_k, pcolor, plabel in pieces:
        # Piece (arc) centre = O_p, with its eccentricity from O_g magnified x60.
        cx = kk * op[0]
        cy = kk * (op[1] - groove[1])
        base = (0.0 if side > 0 else pi) + phi_k
        ts = np.linspace(-half_arc, half_arc, 80)
        arc_x = cx + r_prime * np.cos(base + ts)
        arc_y = cy + r_prime * np.sin(base + ts)
        fx = side * (vane_half + kk * ff * 1e-3)  # flat face at the exaggerated vane gap
        ax.fill(
            list(arc_x) + [fx, fx],
            list(arc_y) + [arc_y[-1], arc_y[0]],  # constant-height flat chord on the vane side
            facecolor=c["bush"],
            ec="#222",
            lw=1.1,
            zorder=4,
        )
        ha = "left" if side > 0 else "right"
        ax.annotate(
            f"곡면 {cf:.2f} µm",
            (cx + r_prime * cos(base), cy + r_prime * sin(base)),
            (side * (rcut + 0.4), 4.4),
            fontsize=8,
            color=c["groove"],
            ha=ha,
            arrowprops={"arrowstyle": "->", "color": c["groove"], "lw": 0.8},
        )
        ax.annotate(
            f"평면 {ff:.2f} µm",
            (fx, cy - 3.0),
            (side * (vane_half + 2.2), -6.4),
            fontsize=8,
            color=c["in"],
            ha=ha,
            arrowprops={"arrowstyle": "->", "color": c["in"], "lw": 0.8},
        )
        ax.plot(cx, cy, "o", color=pcolor, ms=6, mec="white", zorder=6)
        ax.annotate(plabel, (cx, cy), (cx + side * 0.5, cy + 2.4), fontsize=8, color=pcolor, ha=ha)
    ax.plot(ogx, 0, "D", color="#111", ms=6, zorder=6)
    ax.annotate("O_g", (ogx, 0), (ogx - 0.6, 1.5), fontsize=8, ha="right", color="#111")
    ax.annotate("베인", (0, -rcut - 0.2), (0, -rcut - 0.9), fontsize=8.5, color="#eee", ha="center")
    ax.set_xlim(-rcut - 4.5, rcut + 4.5)
    ax.set_ylim(-rcut - 3.4, rcut + 3.2)
    ax.set_xlabel("x (mm, 베인 중선 기준)")
    ax.set_title(
        f"부시 클리어런스 ×{kk} 확대 (개념) — O_p·O_g 치우침", fontsize=10.5, fontweight="bold"
    )


def _build_figure(
    geometry: RotaryGeometry, data: BushOrbitData
) -> tuple[Figure, Callable[[int], None]]:
    """Create the two-panel figure and return ``(fig, update)`` where ``update(j)`` redraws."""

    import matplotlib.pyplot as plt

    _set_korean_font()
    window = _main_window_mm(geometry, data)
    fig = plt.figure(figsize=(12.4, 6.6), dpi=120)
    ax_main = fig.add_axes((0.06, 0.09, 0.50, 0.84))
    ax_inset = fig.add_axes((0.62, 0.09, 0.36, 0.84))

    def update(j: int) -> None:
        draw_main(ax_main, geometry, data, j, window)
        draw_inset(ax_inset, geometry, data, j)

    update(0)
    return fig, update


# ---------------------------------------------------------------------------
# GIF export
# ---------------------------------------------------------------------------


def render_bush_cycle_gif(
    path: str | Path,
    *,
    geometry: RotaryGeometry | None = None,
    data: BushOrbitData | None = None,
    cache_path: str | Path | None = None,
    fps: int = 30,
) -> Path:
    """Render the full-cycle animation to an animated GIF at ``path``."""

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter

    geometry = RotaryGeometry.default() if geometry is None else geometry
    data = compute_orbit_data(geometry, cache_path=cache_path) if data is None else data
    fig, update = _build_figure(geometry, data)
    anim = FuncAnimation(fig, update, frames=data.frames, interval=1000 / fps)  # type: ignore[arg-type]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    anim.save(str(path), writer=PillowWriter(fps=fps))
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Tk viewer (embeds the same figure)
# ---------------------------------------------------------------------------


class BushViewerApp:
    """Tk window embedding the bush figure with a crank-angle slider, play, and Save GIF."""

    def __init__(self, root: tk.Tk, geometry: RotaryGeometry, data: BushOrbitData) -> None:
        import tkinter as tk
        from tkinter import ttk

        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

        self.root = root
        self.data = data
        self.geometry = geometry
        self.frame = 0
        self.playing = False
        root.title("mochi — 부시 클리어런스 / 로터 마우스 애니메이션")

        self.fig, self.update_draw = _build_figure(geometry, data)
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)  # type: ignore[no-untyped-call]
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)  # type: ignore[no-untyped-call]

        bar = ttk.Frame(root, padding=8)
        bar.pack(fill=tk.X)
        self.play_text = tk.StringVar(value="▶ 재생")
        ttk.Button(bar, textvariable=self.play_text, command=self._toggle).pack(side=tk.LEFT)
        self.slider = ttk.Scale(
            bar, from_=0, to=data.frames - 1, orient=tk.HORIZONTAL, command=self._on_slide
        )
        self.slider.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=10)
        self.angle_text = tk.StringVar()
        ttk.Label(bar, textvariable=self.angle_text, width=12).pack(side=tk.LEFT)
        ttk.Button(bar, text="GIF 저장", command=self._save_gif).pack(side=tk.LEFT, padx=(10, 0))

        self._show(0)
        self.root.after(33, self._tick)

    def _show(self, j: int) -> None:
        self.frame = int(j) % self.data.frames
        self.update_draw(self.frame)
        self.angle_text.set(
            f"θ = {_frame_state(self.geometry, self.data, self.frame)['theta_deg']:.0f}°"
        )
        self.canvas.draw_idle()  # type: ignore[no-untyped-call]

    def _on_slide(self, value: str) -> None:
        if not self.playing:
            self._show(int(float(value)))

    def _toggle(self) -> None:
        self.playing = not self.playing
        self.play_text.set("⏸ 일시정지" if self.playing else "▶ 재생")

    def _tick(self) -> None:
        if self.playing:
            self._show(self.frame + 1)
            self.slider.set(self.frame)
        self.root.after(33, self._tick)

    def _save_gif(self) -> None:
        from tkinter import filedialog

        path = filedialog.asksaveasfilename(
            defaultextension=".gif",
            initialfile="bush_clearance_cycle.gif",
            filetypes=[("GIF", "*.gif")],
        )
        if path:
            render_bush_cycle_gif(path, geometry=self.geometry, data=self.data)


def main(argv: list[str] | None = None) -> int:
    """Launch the bush viewer (computes the coupled orbit first, ~minutes)."""

    parser = argparse.ArgumentParser(description="Bush clearance / rotor-mouth animation.")
    parser.add_argument("--gif", type=str, default=None, help="render a GIF to this path and exit")
    parser.add_argument("--cache", type=str, default=None, help="orbit .npz cache (load or write)")
    parser.add_argument(
        "--recompute", action="store_true", help="force re-integration of the orbit"
    )
    args = parser.parse_args(argv)

    geometry = RotaryGeometry.default()
    if args.gif:
        out = render_bush_cycle_gif(args.gif, geometry=geometry, cache_path=args.cache)
        print(f"wrote {out}")
        return 0

    import tkinter as tk

    print("integrating the coupled 9-DOF orbit (a few minutes)...")
    data = compute_orbit_data(geometry, cache_path=args.cache, recompute=args.recompute)
    root = tk.Tk()
    BushViewerApp(root, geometry, data)
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
