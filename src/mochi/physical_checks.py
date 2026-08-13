"""Internal physical-plausibility checks for the film solvers.

These exist because of how the 2026-08 cavitation bug was found. The film pressure field had
been running at a total **absolute** pressure of -11.6 MPa -- tension no lubricant can carry --
for as long as the gas boundary conditions had been in place, and nothing in the code
objected. It surfaced only when an external comparison made someone ask why a friction number
looked large, and the actual evidence, once looked for, was entirely internal: the field
itself was impossible on its own terms.

So the lesson was not that the project needed a better reference. It was that a solver
returning a physically impossible field should say so. The checks here are that voice: cheap,
opt-in assertions over quantities whose valid range follows from physics rather than from
tuning.

They default to **off** in the solvers' hot path -- a coupled orbit evaluates these films tens
of thousands of times -- and are meant to be switched on when a configuration is new or a
result is surprising. :func:`check_pressure_field` is the one that would have caught the bug.
"""

from __future__ import annotations

import numpy as np


class PhysicalityError(ValueError):
    """A solver produced a field or state that cannot exist physically."""


def check_pressure_field(
    pressure: np.ndarray,
    *,
    reference_pressure_pa: float = 0.0,
    vapour_pressure_pa: float = 0.0,
    label: str = "film",
) -> None:
    """Raise when the film pressure implies impossible tension in the lubricant.

    ``pressure`` is the solved field in whatever gauge the caller works in;
    ``reference_pressure_pa`` is the absolute pressure that gauge is referenced to (the
    immersion pressure -- e.g. the 4.0 MPa rotor bore for the bush films, 0 for a film
    referenced to absolute). The absolute field is then ``pressure + reference``, and it may
    not fall below ``vapour_pressure_pa``: a liquid cannot be pulled below its vapour
    pressure, it cavitates instead.

    This is the check that would have caught the 2026-08 bug, where disabling the cavitation
    clamp under a gas boundary let the diverging half-stroke of a reciprocating pad reach
    -11.6 MPa absolute.
    """

    absolute = np.asarray(pressure, dtype=float) + reference_pressure_pa
    worst = float(np.min(absolute))
    if worst < vapour_pressure_pa:
        raise PhysicalityError(
            f"{label}: absolute pressure reaches {worst / 1e6:.3f} MPa, below the vapour "
            f"pressure {vapour_pressure_pa / 1e6:.3f} MPa -- the lubricant would cavitate. "
            f"A cavitation clamp is missing or its floor is set below the vapour pressure."
        )


def check_film_positive(film: np.ndarray, *, label: str = "film") -> None:
    """Raise when the film thickness is non-positive anywhere (surfaces interpenetrating)."""

    film = np.asarray(film, dtype=float)
    worst = float(np.min(film))
    if not np.all(np.isfinite(film)) or worst <= 0.0:
        raise PhysicalityError(
            f"{label}: film thickness reaches {worst * 1e6:.4f} um -- the surfaces "
            f"interpenetrate. The geometry demands contact the model is not carrying."
        )


def check_eccentricity_ratio(
    ratio: float, *, limit: float = 1.0, label: str = "arc film"
) -> None:
    """Raise when an eccentricity ratio exceeds what the clearance geometrically allows.

    A ratio above 1 means the piece centre is further from the groove centre than the
    clearance -- geometrically impossible for rigid surfaces. Reaching it means the load path
    demands a contact reaction; if a hard clamp is absorbing that instead of a contact model,
    the clamp value is deciding the answer (see docs/bush_film_revision_2026-08.md section 5,
    where moving it 0.95 -> 0.999 swings the minimum film by 4900 %).
    """

    if not np.isfinite(ratio) or ratio > limit:
        raise PhysicalityError(
            f"{label}: eccentricity ratio {ratio:.3f} exceeds {limit:.3f} -- the piece is "
            f"further off centre than the clearance allows. A contact model, not a clamp, "
            f"has to carry this."
        )


def check_energy_sign(power_w: float, *, label: str = "friction") -> None:
    """Raise on negative dissipation -- friction cannot return energy to the system."""

    if not np.isfinite(power_w) or power_w < 0.0:
        raise PhysicalityError(
            f"{label}: dissipated power is {power_w:.4f} W. Friction cannot be negative."
        )
