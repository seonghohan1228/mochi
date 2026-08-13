"""Minimal Tecplot ASCII (.dat) writer for exporting raw model data.

The ``results/`` PNGs are for *illustration*; the numbers behind them (film thickness,
per-part forces and moments, validation curves) are exported here as Tecplot ASCII so they
can be post-processed in Tecplot rather than read off a picture. Two zone shapes are
supported, both in ``F=POINT`` layout:

* **1-D point zone** -- one row per sample (e.g. a quantity vs crank angle). Built with
  :func:`point_zone`.
* **2-D ordered zone** -- a structured ``I x J`` grid (e.g. film thickness over
  crank-angle x circumferential position), so Tecplot can draw a contour/surface. Built
  with :func:`ordered_zone`; the point order is I-fastest (Tecplot's ordered convention).

All zones in one file share a single ``VARIABLES`` list (a Tecplot file-level record), so a
file is one dataset with one consistent column set. See ``docs/plotting_conventions.md`` and
PHYSICS.md 4.14.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Zone:
    """One Tecplot zone: a title, POINT-ordered columns, and the ordered dimensions.

    ``columns`` holds one 1-D array per variable, already flattened in Tecplot POINT order
    (I-fastest for a 2-D zone). ``dims`` is ``(I,)`` for a line zone or ``(I, J)`` for a
    structured surface; ``prod(dims)`` must equal every column length.

    ``strand_id`` / ``solution_time`` tag the zone as one frame of a transient series: many
    single-line zones sharing a ``strand_id``, each at its own ``solution_time``, let Tecplot
    step/animate through them (e.g. one film-vs-position line per crank angle).
    """

    title: str
    columns: tuple[np.ndarray, ...]
    dims: tuple[int, ...]
    strand_id: int | None = None
    solution_time: float | None = None


def point_zone(
    title: str, columns: list, *, strand_id: int | None = None, solution_time: float | None = None
) -> Zone:
    """A 1-D point zone from equal-length columns (one per variable).

    Pass ``strand_id`` + ``solution_time`` to make this zone one frame of a transient series.
    """

    cols = tuple(np.asarray(c, dtype=float).ravel() for c in columns)
    if not cols:
        raise ValueError("A zone needs at least one column.")
    n = cols[0].size
    if any(c.size != n for c in cols):
        raise ValueError("All columns in a point zone must have the same length.")
    return Zone(
        title=title, columns=cols, dims=(n,), strand_id=strand_id, solution_time=solution_time
    )


def ordered_zone(
    title: str,
    i_values,
    j_values,
    fields: dict,
    *,
    i_name: str | None = None,
    j_name: str | None = None,
) -> tuple[list[str], Zone]:
    """A 2-D ordered (``I x J``) surface zone; returns ``(variable_names, zone)``.

    ``i_values`` (length I) and ``j_values`` (length J) are the two coordinate axes;
    ``fields`` maps each field name to an ``[I, J]`` array. The returned variable names are
    ``[i_name, j_name, *fields]`` -- pass them (once) as the file's ``variables``. Columns
    are flattened I-fastest so Tecplot reads the structured grid correctly.
    """

    i_values = np.asarray(i_values, dtype=float).ravel()
    j_values = np.asarray(j_values, dtype=float).ravel()
    ni, nj = i_values.size, j_values.size
    ii, jj = np.meshgrid(i_values, j_values, indexing="ij")  # both [I, J]
    columns = [ii.flatten(order="F"), jj.flatten(order="F")]
    names = [i_name or "i", j_name or "j"]
    for name, grid in fields.items():
        arr = np.asarray(grid, dtype=float)
        if arr.shape != (ni, nj):
            raise ValueError(f"Field {name!r} must have shape ({ni}, {nj}), got {arr.shape}.")
        columns.append(arr.flatten(order="F"))
        names.append(name)
    return names, Zone(title=title, columns=tuple(columns), dims=(ni, nj))


def write_dat(
    path: str | Path, variables: list[str], zones: list[Zone], *, title: str | None = None
) -> Path:
    """Write a Tecplot ASCII ``.dat`` file (POINT format) and return its path.

    ``variables`` is the file-level column list (shared by every zone); each zone's column
    count must match ``len(variables)``. Parent folders are created as needed.
    """

    path = Path(path)
    if not zones:
        raise ValueError("Need at least one zone to write.")
    nvar = len(variables)
    for zone in zones:
        if len(zone.columns) != nvar:
            raise ValueError(
                f"Zone {zone.title!r} has {len(zone.columns)} columns but {nvar} variables."
            )
        expected = int(np.prod(zone.dims))
        if any(c.size != expected for c in zone.columns):
            raise ValueError(f"Zone {zone.title!r} column length does not match dims {zone.dims}.")

    path.parent.mkdir(parents=True, exist_ok=True)
    var_record = " ".join(f'"{v}"' for v in variables)
    with path.open("w", encoding="utf-8") as handle:
        if title is not None:
            handle.write(f'TITLE = "{title}"\n')
        handle.write(f"VARIABLES = {var_record}\n")
        for zone in zones:
            if len(zone.dims) == 2:
                header = f'ZONE T="{zone.title}", I={zone.dims[0]}, J={zone.dims[1]}, F=POINT'
            else:
                header = f'ZONE T="{zone.title}", I={zone.dims[0]}, F=POINT'
            if zone.strand_id is not None:
                header += f", STRANDID={zone.strand_id}"
            if zone.solution_time is not None:
                header += f", SOLUTIONTIME={zone.solution_time:.9g}"
            handle.write(header + "\n")
            rows = np.column_stack(zone.columns)
            np.savetxt(handle, rows, fmt="%.9g", delimiter=" ")
    return path
