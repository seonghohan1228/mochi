# mochi

`mochi` is an early-stage Python simulation program for a rotary compressor.
The currently implemented feature is a small prescribed-motion visualization
for checking the cylinder, top-mounted vane, and rotor geometry. It is intended to be
used by the broader simulation; it is not the main numerical solver.

The broader goal remains a simulation whose equations, assumptions, units, and
numerical error can be reviewed independently of its Python implementation.

## Quick start

Requires Python 3.11 or newer.

PowerShell:

```powershell
py -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
python -m pytest
python -m mochi --version
```

Linux or macOS:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
python -m pytest
python -m mochi --version
```

Before opening a pull request, run all local checks:

```bash
ruff format --check .
ruff check .
mypy src/mochi
pytest
```

## Prescribed-motion test GUI

Launch the interactive Tkinter prototype:

```bash
python -m mochi gui
```

The initial inputs are:

| Quantity | Value |
|---|---:|
| Cylinder ID | 77.0 mm |
| Rotor OD | 68.0 mm |
| Eccentricity | 4.5 mm |
| Physical speed | 30 Hz = 1800 rpm |
| Display slow factor | 0.01 |
| Circular cutout radius | 8 mm |
| Cutout center distance | 25 mm from rotor center |
| Vane width | 8 mm |
| Vane-tip calibration | 9 mm from the rotor center when the rotor is at the top |

The inlet marker is 30 degrees to the right of the positive y axis; the
outlet marker is approximately symmetric on the left. These are schematic
port locations, not flow boundary conditions. Angle zero is the top rotor
position and positive angle advances clockwise. The vane length continues to
change while its unfilled outline overlaps the rotor below the circular
cutout. The rotor opening stops at that circle and has tangent fillets at its
two outer lips. The GUI shows prescribed geometry only and does not yet
calculate pressure, contact force, or torque.
The display marks the cylinder center `C`, the moving rotor center `R`, and
the rotor-center orbit as a dashed circle of radius equal to the eccentricity.

## Collaboration model

`main` is the protected, reproducible project history. It should contain only
code that both collaborators have chosen and that passes the agreed checks.
Neither student develops directly on `main`.

When two implementations are developed separately, do not merge both branches
wholesale. Create a short-lived `integration/<scope>` branch from current
`main`, cherry-pick or manually port only the selected commits, verify the
combined solver, and review that integration branch before it reaches `main`.
The unselected implementations remain available in their original branches or
Git history.

See [CONTRIBUTING.md](CONTRIBUTING.md) for commands, branch names, conflict
handling, and the exact selection workflow.

## Documentation map

- [PHYSICS.md](PHYSICS.md): governing equations, coordinates, units,
  assumptions, numerical method, and validation contract.
- [ARCHITECTURE.md](ARCHITECTURE.md): module boundaries and data flow.
- [CONTRIBUTING.md](CONTRIBUTING.md): the two-person Git and review workflow.
- [PLAN.md](PLAN.md): staged development plan and unresolved decisions.
- [CHANGES.md](CHANGES.md): user-visible changes and why they were made.
- [AGENTS.md](AGENTS.md): repository instructions for coding agents.

If code and documentation disagree, treat that as a defect. Physics changes
must update `PHYSICS.md`; workflow or public-interface changes must update the
corresponding guide in the same pull request.

## Project layout

```text
src/mochi/             installable Python package
  cli.py               command-line entry point
  kinematics.py        prescribed rotary mechanism geometry and motion
  gui.py               Tkinter test GUI and animation
  __main__.py          `python -m mochi` entry point
tests/                 fast unit and regression tests
.github/workflows/     checks run for branches and pull requests
AGENTS.md               instructions for automated coding agents
PHYSICS.md              model and numerical-method contract
ARCHITECTURE.md         intended package boundaries
CONTRIBUTING.md         collaboration and Git workflow
PLAN.md                 development stages and open decisions
CHANGES.md              append-only change history
```

Generated simulations belong under `results/` and are ignored by Git. Small,
reviewed reference data may later be committed under `tests/data/` with its
origin, units, and generation procedure documented.

## Current scope

The visualization prescribes the rotor-center orbit and keeps the internal
circular-cutout center on the stationary top-vane axis. It does not calculate
pressure, force, contact, or torque. Those models belong to the broader
simulation and remain outside the accepted implementation until their equations
and validation cases are documented in `PHYSICS.md`.

## License

This project is released under the terms in [LICENSE](LICENSE).
