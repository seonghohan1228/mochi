# Figure / plotting conventions

Rules for every image under `results/`, produced by
`scripts/generate_results.py`. These are standing **user requirements** — keep
them in every session. (See also the one-line reminder in AGENTS.md → *`results/`
figures*.)

## 1. One plot per figure

Each figure file is a **single plot (a single Matplotlib axes)** — never a
multi-panel montage (`plt.subplots(1, 3, ...)`, side-by-side panels A/B/C, etc.).
If a topic needs several views, produce **several files**, one plot each. The
`FIGURES` registry comment already states "one single-axes figure per file"; do
not regress to montages.

## 2. Group same-topic figures in a folder

Related figures share one `results/<topic>/` folder (e.g. `bush_film/`,
`bearing_load/`, `assembly/`). A new topic gets a new folder; the filename names
the specific view within the topic.

## 3. Annotations must not obscure the graph

- Legends go **outside** the axes (use the `_legend_above` helper) or in a
  genuinely empty corner — never over the data.
- Value labels / callouts sit in empty regions with a short leader arrow to the
  point; they must not cover curves, markers, or each other.
- Reserve margin: widen `xlim`/`ylim` or place text boxes outside the data area
  so nothing overlaps. If a label has nowhere clear to go, the figure is trying
  to show too much — split it (rule 1).

## 4. Grid / convergence tests → grid-count vs relative-error graph

Any grid- or mesh-convergence study is presented as a **grid count (x) vs relative
error (y)** graph, with:

- **one figure per quantity/channel** — never several channels/films on one axes;
  split them (e.g. a separate convergence figure for each oil film / friction
  channel);
- **monochrome (black)** line and markers — the message is the convergence *trend*,
  not a categorical colour;
- error against the finest grid (`|value - value_finest| / |value_finest|`), log
  y-axis (log-log when the grid spans decades);
- report the converged value and the grid where the relative error drops below the
  target (e.g. 1 %).

Group the figures in a `results/convergence/` folder and keep the numeric table in a
summary file (`docs/*_convergence.md` or a `.txt`).

## Rationale

The gallery is read one image at a time. Multi-panel montages and annotations
that sit on top of the curves have been flagged repeatedly. Figures are
regenerated from `generate_results.py`, never hand-edited, so these rules live in
the renderers.
