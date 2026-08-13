"""Tests for the Tecplot ASCII (.dat) writer (:mod:`mochi.tecplot`)."""

import numpy as np
import pytest

from mochi.tecplot import ordered_zone, point_zone, write_dat


def _read(path):
    return path.read_text(encoding="utf-8").splitlines()


def test_point_zone_roundtrips_columns(tmp_path):
    theta = np.linspace(0.0, 360.0, 5)
    h = np.arange(5.0)
    path = write_dat(
        tmp_path / "line.dat", ["theta_deg", "h_um"], [point_zone("s", [theta, h])], title="t"
    )
    lines = _read(path)
    assert lines[0] == 'TITLE = "t"'
    assert lines[1] == 'VARIABLES = "theta_deg" "h_um"'
    assert lines[2] == 'ZONE T="s", I=5, F=POINT'
    body = np.array([[float(x) for x in row.split()] for row in lines[3:] if row.strip()])
    assert body.shape == (5, 2)
    np.testing.assert_allclose(body[:, 0], theta)
    np.testing.assert_allclose(body[:, 1], h)


def test_ordered_zone_is_i_fastest(tmp_path):
    # A 3x2 field flattened I-fastest: for each j, i runs 0..2.
    field = np.arange(6.0).reshape(3, 2)  # [[0,1],[2,3],[4,5]]
    names, zone = ordered_zone(
        "m", [0.0, 90.0, 180.0], [0.0, 10.0], {"h": field}, i_name="phi", j_name="theta"
    )
    assert names == ["phi", "theta", "h"]
    path = write_dat(tmp_path / "map.dat", names, [zone])
    lines = _read(path)
    assert lines[0] == 'VARIABLES = "phi" "theta" "h"'  # no TITLE line when omitted
    assert lines[1] == 'ZONE T="m", I=3, J=2, F=POINT'
    body = np.array([[float(x) for x in row.split()] for row in lines[2:] if row.strip()])
    # First block is theta=0 with phi = 0,90,180 and h = column 0 of the field (0,2,4).
    np.testing.assert_allclose(body[:3, 0], [0.0, 90.0, 180.0])
    np.testing.assert_allclose(body[:3, 1], [0.0, 0.0, 0.0])
    np.testing.assert_allclose(body[:3, 2], [0.0, 2.0, 4.0])
    # Second block is theta=10 with h = column 1 (1,3,5).
    np.testing.assert_allclose(body[3:, 2], [1.0, 3.0, 5.0])


def test_point_zone_rejects_ragged_columns():
    with pytest.raises(ValueError):
        point_zone("bad", [np.arange(3.0), np.arange(4.0)])


def test_ordered_zone_rejects_wrong_field_shape():
    with pytest.raises(ValueError):
        ordered_zone("bad", [0.0, 1.0, 2.0], [0.0, 1.0], {"h": np.zeros((2, 2))})


def test_transient_strand_headers(tmp_path):
    # Two frames of a strand: each a 1-D line at its own solution time, shared strand id.
    zones = [
        point_zone("f0", [np.arange(3.0), np.arange(3.0)], strand_id=1, solution_time=0.0),
        point_zone("f1", [np.arange(3.0), np.arange(3.0) + 1], strand_id=1, solution_time=0.5),
    ]
    path = write_dat(tmp_path / "strand.dat", ["x", "h"], zones)
    headers = [line for line in _read(path) if line.startswith("ZONE")]
    assert len(headers) == 2
    assert 'ZONE T="f0", I=3, F=POINT, STRANDID=1, SOLUTIONTIME=0' == headers[0]
    assert "STRANDID=1, SOLUTIONTIME=0.5" in headers[1]


def test_write_dat_checks_variable_count(tmp_path):
    zone = point_zone("s", [np.arange(3.0), np.arange(3.0)])
    with pytest.raises(ValueError):
        write_dat(tmp_path / "x.dat", ["only_one_name"], [zone])
