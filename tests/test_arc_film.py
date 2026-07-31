"""Partial-arc axial-uniform 1-D film force (curved swing-bush piece film)."""

import math

import pytest

from mochi.arc_film import arc_film_force
from mochi.journal_bearing import JOURNAL_CLEARANCE_M, JOURNAL_LENGTH_M, JOURNAL_RADIUS_M
from mochi.slider_film import flat_slider_film

R, L, C = JOURNAL_RADIUS_M, JOURNAL_LENGTH_M, JOURNAL_CLEARANCE_M


def _arc(ex, ey, edx, edy, omega, half=math.pi, center=0.0, n=2001):
    return arc_film_force(
        ex,
        ey,
        edx,
        edy,
        omega,
        arc_center_rad=center,
        arc_half_span_rad=half,
        radius_m=R,
        length_m=L,
        clearance_m=C,
        n_beta=n,
    )


def test_narrow_arc_matches_the_flat_slider() -> None:
    """A shallow arc is a flat wedge: its load matches the analytic slider film.

    For ``e_x = 0`` a small half-span ``D`` about ``beta = 0`` gives a nearly straight
    wedge ``h ~ c - (e_y/r) x`` with mean surface speed ``Omega r``; the arc load must
    therefore agree with :func:`flat_slider_film` fed the matched wedge
    (``gamma = -e_y/r``, ``U = 2 Omega r``, length ``2 r D``).
    """

    half, e_y, omega = 0.12, 0.3 * C, 95.0
    arc = _arc(0.0, e_y, 0.0, 0.0, omega, half=half)
    slider = flat_slider_film(
        0.0,
        -e_y / R,
        0.0,
        0.0,
        2.0 * omega * R,
        length_m=2.0 * R * half,
        height_m=L,
        clearance_m=C,
        n_s=4001,
    )
    assert abs(arc.force_x_n) == pytest.approx(slider.normal_force_n, rel=3.0e-2)
    assert abs(arc.force_y_n) < 0.1 * abs(arc.force_x_n)  # tangential is second order


def test_force_is_restoring_along_the_line_of_centres() -> None:
    """Eccentricity toward +x gives a restoring force with a -x component."""

    f = _arc(0.5 * C, 0.0, 0.0, 0.0, 95.0)
    assert f.force_x_n < 0.0  # pushes the piece back toward centre
    assert f.min_film_thickness_m == pytest.approx(0.5 * C, rel=1e-9)


def test_partial_arc_carries_less_than_full_arc() -> None:
    """A half arc supports a smaller load than the full ring at the same eccentricity."""

    full = _arc(0.5 * C, 0.0, 0.0, 0.0, 95.0, half=math.pi)
    part = _arc(0.5 * C, 0.0, 0.0, 0.0, 95.0, half=math.pi / 2)
    assert math.hypot(part.force_x_n, part.force_y_n) < math.hypot(full.force_x_n, full.force_y_n)


def test_squeeze_generates_force_at_zero_speed() -> None:
    """A closing film (piece pushed into the arc, edot_x > 0) carries load with no speed."""

    f = _arc(0.3 * C, 0.0, 0.01, 0.0, 0.0, half=math.pi / 2)
    assert math.hypot(f.force_x_n, f.force_y_n) > 0.0


def test_tangential_eccentricity_gives_tangential_force_on_partial_arc() -> None:
    """On a partial arc a tangential (y) offset produces a y-force (broken symmetry)."""

    f = _arc(0.0, 0.4 * C, 0.0, 0.0, 95.0, half=math.pi / 3)
    assert abs(f.force_y_n) > 1.0e-6


def test_force_is_grid_converged() -> None:
    """The load is stable under grid refinement (the line solve has converged)."""

    coarse = _arc(0.4 * C, 0.1 * C, 0.0, 0.0, 95.0, half=math.pi / 2, n=401)
    fine = _arc(0.4 * C, 0.1 * C, 0.0, 0.0, 95.0, half=math.pi / 2, n=3001)
    assert coarse.force_x_n == pytest.approx(fine.force_x_n, rel=2.0e-3)
    assert coarse.force_y_n == pytest.approx(fine.force_y_n, rel=2.0e-3)


def test_contact_raises() -> None:
    with pytest.raises(ValueError):
        _arc(C, 0.0, 0.0, 0.0, 95.0)  # eccentricity == clearance -> zero film


def test_invalid_inputs_raise() -> None:
    with pytest.raises(ValueError):
        arc_film_force(
            0.0,
            0.0,
            0.0,
            0.0,
            95.0,
            arc_center_rad=0.0,
            arc_half_span_rad=0.0,
            radius_m=R,
            length_m=L,
            clearance_m=C,
        )
    with pytest.raises(ValueError):
        arc_film_force(
            0.0,
            0.0,
            0.0,
            0.0,
            95.0,
            arc_center_rad=0.0,
            arc_half_span_rad=math.pi,
            radius_m=R,
            length_m=L,
            clearance_m=C,
            n_beta=4,
        )
