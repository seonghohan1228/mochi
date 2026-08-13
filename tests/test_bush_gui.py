"""Smoke tests for the bush clearance / rotor-mouth animation module.

These exercise the drawing/GIF path on a tiny *synthetic* orbit so the expensive coupled
9-DOF integration is never run in the test suite. Skipped when matplotlib is absent.
"""

import numpy as np
import pytest

pytest.importorskip("matplotlib")
import matplotlib

matplotlib.use("Agg")

from mochi import bush_gui as bg
from mochi.kinematics import RotaryGeometry


def _synthetic(n: int = 3) -> bg.BushOrbitData:
    """A finite, physically plausible few-frame orbit (positions in metres, films in metres)."""

    z = np.zeros(n)
    return bg.BushOrbitData(
        crank=np.linspace(0.5, 2.5, n),
        ejx=z + 1e-6,
        ejy=z + 1e-6,
        rdev=z + 1e-4,
        inx=z + 1e-4,
        iny=z + 0.028,
        inphi=z,
        outx=z - 1e-4,
        outy=z + 0.028,
        outphi=z,
        inc=z + 0.5e-6,
        inf=z + 2e-6,
        outc=z + 20e-6,
        outf=z + 8e-6,
        cj=15e-6,
    )


def test_orbit_data_reports_frame_count() -> None:
    assert _synthetic(4).frames == 4


def test_frame_state_is_finite_and_wrapped() -> None:
    state = bg._frame_state(RotaryGeometry.default(), _synthetic(), 0)
    assert np.isfinite(state["groove"]).all()
    assert 0.0 <= state["theta_deg"] < 360.0
    assert np.isfinite([state["in_c"], state["out_c"], state["in_f"], state["out_f"]]).all()


def test_figure_builds_and_updates() -> None:
    import matplotlib.pyplot as plt

    data = _synthetic()
    fig, update = bg._build_figure(RotaryGeometry.default(), data)
    for j in range(data.frames):
        update(j)  # must not raise for any frame
    plt.close(fig)


def test_gif_renders_to_file(tmp_path) -> None:
    out = bg.render_bush_cycle_gif(
        tmp_path / "cycle.gif", geometry=RotaryGeometry.default(), data=_synthetic(), fps=5
    )
    assert out.exists() and out.stat().st_size > 0
