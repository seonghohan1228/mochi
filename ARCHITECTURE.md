# Architecture

`mochi` should keep compressor physics independent from integration algorithms,
configuration formats, and output. The current code includes prescribed
kinematics and a Tkinter test GUI; the remaining structure below is the
intended direction, not a claim that every module already exists.

## Data flow

```text
case input
    -> validated configuration and immutable model data
    -> initial state
    -> simulation orchestrator
         -> kinematics and force models
         -> integrator through derivative/residual interface
    -> result object with diagnostics
    -> optional file output and plots
```

Dependencies point inward toward state and model interfaces. Physics does not
read files or parse command-line arguments. Numerical integrators do not import
rotary-compressor force implementations.

## Intended package boundaries

Add modules only when their behavior is implemented and tested:

```text
src/mochi/
  __init__.py          deliberately small public API
  cli.py               command parsing; no physics equations
  kinematics.py        current prescribed rotary mechanism
  gui.py               current Tkinter visualization and controls
  config.py            validated case configuration and unit conversion
  state.py             state and result data structures
  simulation.py        assembles a case and coordinates a solve
  physics/
    forces.py          force protocol and total-force assembly
    contact.py         future contact and constraint laws
    ...                one focused module per accepted physical model
  numerics/
    integrators.py     solver adapters and convergence/failure diagnostics
  io/
    cases.py           versioned case serialization
    results.py         versioned result serialization
```

Tests mirror the package boundaries. Small reference data belongs in
`tests/data/` only with provenance and units.

## Core interfaces

The solver should converge on a few explicit concepts:

- **Model data:** validated, immutable geometry and material parameters.
- **State:** time, generalized position, and generalized velocity with known
  array shapes and SI units.
- **Force contribution:** a pure evaluation from time, state, and model data to
  a force on the modeled body in the global frame.
- **Dynamics function:** assembles forces and returns a derivative or residual.
- **Integrator:** advances a generic dynamics function and reports numerical
  diagnostics without knowing compressor details.
- **Result:** state history, individual force histories, metadata, and terminal
  status.

Prefer dataclasses and NumPy arrays at these boundaries. Validate shapes and
finite values at public entry points. Do not pass unstructured dictionaries
through the numerical core.

## Configuration and units

The numerical core uses SI units and radians. Unit conversion belongs at input
and output boundaries. Configuration must be validated once before integration;
the time loop must not repeatedly parse strings or infer missing units.

A serialized case format should receive an explicit schema version before it
is treated as stable. Defaults must be visible and physically defensible. A
missing model is different from a model configured with a zero coefficient.

## Results and reproducibility

Generated output belongs under `results/` and is not committed. Each run should
eventually record the input case, Git revision, package version, active models,
solver settings, status, and output variables. Never overwrite a trusted
reference case without explicit review.

## Error handling

Invalid input raises a specific error before the solve. A numerical failure
returns or raises a result that preserves the reason and last valid state. NaN,
infinite values, contact penetration beyond an accepted tolerance, and solver
nonconvergence must not pass silently.

## Performance policy

Begin with clear vectorized NumPy/SciPy code and measure representative cases
before optimizing. Keep a correct reference implementation when adding a
faster algorithm. Parallelism, compilation, or accelerators require a measured
bottleneck and numerical-equivalence tests.

## Portability

Windows and Linux are first-class development targets. Use `pathlib.Path` for
paths, avoid platform-specific shell assumptions in Python, and keep output
paths relative to an explicit case or project root rather than the process's
accidental working directory.
