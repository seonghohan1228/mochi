# Contributing

This project uses a review-first, two-person workflow. The purpose is to let
both collaborators explore independently while keeping the accepted solver
small, understandable, and reproducible.

## Terms used here

- **`main` branch:** the protected official project history.
- **working branch:** one student's feature, fix, or experiment.
- **accepted solver:** the selected combination being evaluated. Before
  approval it lives on an integration branch; after approval it lives on
  `main`.
- **integration branch:** a short-lived branch used to assemble selected
  changes. It is not a permanent `develop` branch.

This distinction avoids treating ?merge into the main code? as ?merge every
development branch into `main`.? Selection happens first; `main` receives only
the verified result.

## Protect `main`

Recommended repository settings:

- Require a pull request before merging.
- Require the other collaborator's approval.
- Require the Linux and Windows CI checks.
- Require conversations to be resolved.
- Block force pushes and branch deletion.
- Prefer squash merging for a noisy exploratory branch and fast-forward or
  regular merging for an already curated integration branch.

Direct commits to `main` are reserved for a genuine repository recovery, not
routine convenience. `main` should always be installable, testable, and tied
to documented physics.

## Branch names

Use lowercase, short, descriptive names:

```text
feature/<owner>/<scope>       feature/sh/kinematics
experiment/<owner>/<scope>    experiment/yk/contact-penalty
fix/<owner>/<scope>           fix/sh/angle-sign
integration/<scope>           integration/contact-model
docs/<owner>/<scope>          docs/yk/force-conventions
```

The owner segment can be initials or a stable short username. Avoid permanent
personal branches; create a fresh branch for each question or change.

## Start and develop a task

```bash
git switch main
git pull --ff-only
git switch -c feature/<owner>/<scope>
```

During development:

1. Change one coherent idea at a time.
2. Add or update its tests and documentation in the same branch.
3. Inspect the staged diff before each commit.
4. Commit a logical selection unit with an imperative message.

```bash
git add <intentional-files>
git diff --staged
git commit -m "Add planar bushing state"
```

Good commits make later selection easy. Do not combine a force model, unrelated
formatting, and generated results in one commit.

Before review, update a private branch from `main`:

```bash
git fetch origin
git rebase origin/main
```

Rebase only a private branch. If another person is already using the branch,
merge `origin/main` into it instead and do not rewrite their history.

## Select changes from parallel work

Suppose both students produced candidate implementations. Review them against
the same acceptance cases and then:

```bash
git switch main
git pull --ff-only
git switch -c integration/<scope>
git cherry-pick <selected-commit-1>
git cherry-pick <selected-commit-2>
```

Use cherry-pick when a whole commit is accepted. If only part of a commit is
accepted, manually port the chosen lines on the integration branch, test them,
and create a new commit. In that commit or pull request, record the source
branch and commit so authorship and reasoning remain traceable.

Do not merge both candidate branches and then try to remove rejected work. That
usually retains accidental coupling and makes review harder. If the two
approaches edit the same equations, appoint one integrator for that task; the
other collaborator reviews the combined result.

On the integration branch:

1. Resolve conflicts according to `PHYSICS.md`, not merely by choosing the code
   that runs.
2. Run analytic/reference cases for each candidate and the combined solver.
3. Run formatting, linting, type checking, and the full test suite.
4. Update `PHYSICS.md`, `PLAN.md`, and `CHANGES.md` as appropriate.
5. Open one pull request from the integration branch to `main` and list exactly
   which commits or ideas were selected and rejected.

## Conflicts

Never resolve a physics conflict by accepting ?ours? or ?theirs? without
checking signs, units, frames, and assumptions. For a nontrivial conflict:

1. Write down the intended equation or interface.
2. Add a test that distinguishes the alternatives.
3. Resolve the code to satisfy the agreed contract.
4. Review the final diff, not only the conflict markers.

Abort the operation and ask the other collaborator if unrelated work appears
or the intended model is unclear.

## Pull request checklist

- The branch has one clear purpose.
- Selected source branches/commits are identified.
- Equations, frames, units, and assumptions are documented.
- New behavior has a focused test and appropriate validation evidence.
- `ruff format --check .`, `ruff check .`, `mypy src/mochi`, and `pytest` pass.
- Linux and Windows CI pass.
- No generated results, virtual environments, secrets, or unrelated edits are
  included.
- User-facing and physics documentation match the implementation.

After merge, both collaborators update local `main`. Delete short-lived feature
and integration branches once their useful comparison history is no longer
needed. Use tags such as `v0.1.0` only for states worth reproducing or sharing.
