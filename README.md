# mochi

`mochi` is an early-stage Python numerical solver for the motion of a bushing
in a scroll compressor. The repository currently provides the project shell,
collaboration rules, and model contract; the compressor geometry, force laws,
and production time integrator still need to be agreed and validated.

The goal is a solver whose equations, assumptions, units, and numerical error
can be reviewed independently of its Python implementation.

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
src/mochi/             installable Python package and command-line entry point
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

The initial solver target is planar rigid-body translation of the bushing from
a defined initial state under explicitly assembled forces. Rotation, contact,
friction, hydrodynamic-film forces, thermal effects, and compressor-cycle
coupling are not assumed to be part of the first model until their equations
and validation cases are accepted in `PHYSICS.md`.

## License

This project is released under the terms in [LICENSE](LICENSE).
