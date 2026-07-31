"""Indicated work and power from the port-timed pressure cycle.

The indicated work of a chamber is the area of its closed pressure-volume loop
over one revolution, ``W = -oint p dV``; the indicated power is that work times
the shaft frequency, ``P = W f``. Both use the accepted models directly: the
pressure comes from the port-timed true-geometry rule of PHYSICS.md section 3.4
(:func:`mochi.chambers.port_timed_pressures`) and the volume from the cached
true-geometry volume trace (:class:`mochi.chambers.CycleTrace`).

Sign convention: work done *on* the gas by the mechanism is positive, so the
loop integral carries a leading minus sign and a compressor returns a positive
indicated work and power. See PHYSICS.md section 3.5.

This is a derived performance quantity built on the prescribed pressure and
volume boundary conditions; it is not a transient thermodynamic solution and it
neglects leakage, heat transfer, and mechanical losses (the last are estimated
separately by :mod:`mochi.bush_film`).
"""

from dataclasses import dataclass
from math import pi

from mochi.chambers import (
    CycleTrace,
    build_cycle_trace,
    port_timed_pressures,
    seal_over_half_angle_rad,
)
from mochi.kinematics import RotaryGeometry

FULL_TURN_RAD = 2.0 * pi


@dataclass(frozen=True, slots=True)
class IndicatedWork:
    """Indicated work per revolution and power for both chambers (SI units)."""

    compression_work_j: float
    suction_work_j: float
    net_work_j: float
    power_w: float
    frequency_hz: float
    samples: int


def indicated_work_j(
    geometry: RotaryGeometry,
    *,
    samples: int = 720,
    trace: CycleTrace | None = None,
) -> IndicatedWork:
    """Return the indicated work and power over one revolution.

    The pressure-volume loop is integrated with the trapezoidal rule on
    ``samples`` equally spaced crank angles. ``W = -oint p dV`` is evaluated for
    the compression (discharge) chamber and the suction chamber; ``net_work_j``
    is their sum and ``power_w`` multiplies it by ``geometry.frequency_hz``.
    Positive means net work input to the gas (compressors return positive).

    Integration stops at the seal-over entry angle (``2*pi`` minus the seal-over
    half-angle) rather than the full turn: inside the seal-over window the two
    chambers merge and the true-geometry split becomes degenerate, so the
    volume trace balloons there. The compression branch (max to clearance
    volume) and the suction branch (clearance to max volume) over ``[0,
    seal-over entry]`` already form the closed indicator loop; the excluded
    clearance-gas re-expansion across the merge window carries negligible work.
    """

    geometry.validate()
    if samples < 8:
        raise ValueError("Indicated work needs at least eight crank-angle samples.")

    trace = build_cycle_trace(geometry) if trace is None else trace

    seal_over_entry = FULL_TURN_RAD - seal_over_half_angle_rad(geometry)
    angles = tuple(seal_over_entry * index / samples for index in range(samples + 1))
    pressures = tuple(port_timed_pressures(geometry, angle, trace=trace) for angle in angles)
    compression_p = tuple(p.compression_pressure_pa for p in pressures)
    suction_p = tuple(p.suction_pressure_pa for p in pressures)
    compression_v = tuple(trace.compression_at(angle) for angle in angles)
    suction_v = tuple(trace.suction_at(angle) for angle in angles)

    compression_work = 0.0
    suction_work = 0.0
    for index in range(samples):
        compression_work -= (
            0.5
            * (compression_p[index] + compression_p[index + 1])
            * (compression_v[index + 1] - compression_v[index])
        )
        suction_work -= (
            0.5
            * (suction_p[index] + suction_p[index + 1])
            * (suction_v[index + 1] - suction_v[index])
        )

    net_work = compression_work + suction_work
    return IndicatedWork(
        compression_work_j=compression_work,
        suction_work_j=suction_work,
        net_work_j=net_work,
        power_w=net_work * geometry.frequency_hz,
        frequency_hz=geometry.frequency_hz,
        samples=samples,
    )
