# Grid convergence (improved model): friction channels and bush clearance

> **Which channels have a physical spatial grid.** Only the **two bush films** are solved on a
> spatial mesh: the flat slider (`n_s`) and the curved arc (`n_beta`), each a 1-D Reynolds line
> discretisation. **A genuine grid-convergence study exists only for these two.** The **journal**
> force is the **Ocvirk closed form** (`short_bearing_force`) and the **seal** is an analytical
> **Hertz line contact** (`hertz_line_contact_force_n`) — both evaluated *exactly* at each crank
> angle, with **no spatial mesh to refine**. Their only numerical grid is the **cycle time sampling**
> (`samples`), which is a *time condition*, not a physical grid. The `samples`/`grid_samples` sweeps
> below for the journal and seal are therefore **time-sampling checks, not spatial convergence** —
> labelled as such. (`grid_samples` is likewise the gas-BC prescription grid over the cycle, a time
> grid.) The physically meaningful convergence lives in the two bush-film sections and the bush
> **clearance** section.

## Friction channels

Grid/mesh convergence of the four friction channels vs **grid count**, with the
time-integration length **fixed at `revolutions=8`** (warm-start, where every channel has
reached its steady periodic orbit — see the steady-state note). One grid parameter is varied
at a time, in a **doubling (×2) sequence with at least four grid levels** per channel. Model:
the **improved** bush orbit (gas-pressure film BC + flat-film EHL/mixed + curved-film shear
moment, seal on); journal from `integrate_rotor_orbit`. Relative error is `|v - v_ref| / |v_ref|`
with `v_ref` the **Richardson-extrapolated (h→0) limit** from the three finest levels
(observed-order, grid ratio 2) — so **every grid level, including the finest, has a non-zero
error and appears on the log axis** (against the finest grid the finest point would be 0 and
vanish). For the non-converging curved channel `v_ref` is the sample mean (a noise floor).
Figures (one black plot per channel, log y, every point value-labelled):
`results/convergence/{journal, bush_curved, seal, bush_flat}_friction.png`.
Reference: `samples=180, grid_samples=180, n_beta=121, n_s=81`, `revolutions=8`,
`FLAT_ASPERITY_PRESSURE_PA=50 MPa`.

## Summary

| channel | grid (×2, ≥4 levels) | value (W) | grid-convergence |
|---|---|---|---|
| rotor–crank (journal) | samples 90/180/360/720 | 10.18 | rel err → <0.05 % |
| rotor–cylinder (seal) | samples 90/180/360/720 | 11.0 | rel err → <0.4 % |
| bush–vane (flat) | n_s 41/81/161/321 | 53.1 | rel err → <0.05 % |
| bush–rotor (curved) | grid_samples 180…1440 | ~0.1–0.19 | **does NOT converge — ±28 % noise** (negligible, <0.3 %) |

Three channels converge monotonically below 1 %. The **curved** channel does **not** converge on
any grid: it depends on the small velocity difference `φ̇_rotor − φ̇_piece` (whose two states carry
a deliberately loose `atol=1e-6`, see the steady-state note), so the tiny `μu²/h` loss scatters
±28 % about ~0.14 W with the gas-BC grid and no monotone trend. It is ~0.3 % of the total budget,
so this is an honest un-pinnable floor, not a modelling failure — and it does not perturb the other
three channels or the friction budget.

## Journal (rotor–crank), W — samples 90/180/360/720

| samples | 90 | 180 | 360 | 720 |
|---|---|---|---|---|
| value | 10.153 | 10.174 | 10.184 | 10.189 |
| rel err vs limit (10.194 W, p≈1.0) | 0.40 % | 0.20 % | 0.098 % | 0.049 % |

Monotone, error halving with each doubling; below 0.2 % by `samples≥180`. Gas-grid-independent.

## Bush–rotor (curved film), W — does NOT converge

The curved figure plots the value against **`grid_samples`** (the gas-BC prescription grid), which
is the only knob that moves it — held at `n_beta=241, samples=360, rev=8`:

| grid_samples | 180 | 360 | 720 | 1440 |
|---|---|---|---|---|
| curved (W) | 0.1069 | 0.1848 | 0.1254 | 0.1453 |

Refining the gas grid **8×** leaves the value bouncing ±28 % about ~0.14 W with **no monotone
trend** — it does not converge. Two other grids were ruled out as the cause:

- **`samples` (cycle-average points): no effect.** `(samples 180→360, grid_samples 180)` gives
  0.1074 → 0.1069 — unchanged. The value is a post-integration `np.mean`, not a time-resolution issue.
- **`n_beta` (curved line grid): no effect.** At `grid_samples=180` the n_beta sweep 61/121/241/481
  gives 0.1133/0.1068/0.1074/0.1140 (±3 %) — the ±3 % there was just the small residual at one
  under-resolved gas grid, not a convergence trend.

**Root cause.** The curved friction is `μ u² / h` with `u = |φ̇_rotor − φ̇_piece|·r` — a *difference*
of two of the nine velocity states. Those velocities carry a deliberately **loose `atol=1e-6`** (so
BDF can step over the stiff squeeze mode; positions stay tight at 1e-10/1e-11). Every other channel
depends only on positions/films and is clean; the curved loss is the one quantity built from a small
velocity difference, so a sub-percent gas-grid perturbation of the orbit swings this ~0.14 W signal
±28 %. Pinning it would need a tight velocity `atol`, which re-stiffens the integrator — not worth it
for a channel that is **~0.3 % of the total friction**. Reported honestly as an un-pinnable floor.

## Rotor–cylinder (seal), W — samples 90/180/360/720

| samples | 90 | 180 | 360 | 720 |
|---|---|---|---|---|
| value | 11.261 | 11.135 | 11.066 | 11.027 |
| rel err vs limit (10.976 W, p≈0.8) | 2.6 % | 1.4 % | 0.82 % | 0.46 % |

Monotone; drops below 1 % between `samples=180` and 360. (Seal is rev-independent — the 720 level
was taken at rev=2 for speed, consistent with the rev-8 trend.)

## Bush–vane (flat film), W — n_s 41/81/161/321 (rev=8, warm-start)

| n_s | 41 | 81 | 161 | 321 |
|---|---|---|---|---|
| value | 53.145 | 53.130 | 53.123 | 53.121 |
| rel err vs limit (53.120 W) | 0.047 % | 0.018 % | 0.0053 % | 0.0015 % |

Monotone and already **<0.05 % at the coarsest slider grid** (`n_s=41`), confirming the flat-film
loss is grid-insensitive once warm-started. Crucially the flat loss is **insensitive to the same
`grid_samples` that wrecks the curved channel**: at `n_beta=241, samples=360`, `grid_samples`
180/360/720/1440 gives 53.256/53.422/53.438/53.455 W — a **monotone +0.4 %** drift, converged. The
53 W bush loss and the friction budget are therefore safe; only the tiny velocity-difference curved
term is noise-dominated.

## Bush curved CLEARANCE (position) — physical grid `n_beta` 61/121/241/481

Because *predicting the bush position* matters more than its friction, the curved-film gap (the
quantity that fixes where each piece sits) is grid-converged on its **physical** mesh `n_beta`, with
time conditions unified (`samples=360, grid_samples=360, rev=8`). Per-cycle max/min per piece, in µm:

| n_beta | IN min | IN max | OUT min | OUT max |
|---|---|---|---|---|
| 61 | 0.300\* | 1.876 | 16.393 | 32.207 |
| 121 | 0.300\* | 1.818 | 16.180 | 32.180 |
| 241 | 0.300\* | 1.845 | 16.399 | 32.209 |
| 481 | 0.339\* | 1.932 | 16.212 | 32.180 |
| grid scatter (vs mean) | *clamp | ±3 % | ±0.7 % | **±0.05 %** |

Figures (physical grid, log rel-error vs the mean): `results/convergence/clr_curved_{out_max, out_min,
in_max, in_min}.png`.

- **Position is grid-robust** — unlike the ±28 % curved *friction*, the curved *clearance* barely
  moves with `n_beta`: OUT max ±0.05 %, OUT min ±0.7 %, IN max ±3 % (the thin-film piece). The
  clearance rides the *positions* (tight `atol` 1e-10/1e-11), not the loose velocity difference, so it
  converges where the friction cannot.
- **\*IN min is a clamp, not a physical film.** The IN (suction-side) piece is driven to the
  eccentricity clamp `_MAX_ARC_ECC_RATIO=0.99`, so its min gap sits at the floor
  `curved_gap·(1−0.99) = 30 µm × 0.01 = 0.30 µm` (near metal contact). This is a **known model limit**
  (§4.11): the true IN min gap needs a vane-bush **EHL/contact** treatment; the reported 0.30 µm is the
  clamp, held flat across `n_beta` until `n_beta=481` lifts it to 0.339 µm. The OUT piece (discharge)
  keeps a healthy 16–32 µm film and converges cleanly.

## Steady-state note (revolutions is fixed, not a grid)

`revolutions` is the time-integration length (a steady-state parameter), not a mesh, so it is held
fixed at 8 for the grid study. Three channels reach steady state within a few revolutions. The
**flat film** is the exception on a **cold start** (film at the 10 um reference): it drains slowly
toward its ~1.3 um contact equilibrium, so the loss climbs 21.6 (rev 2) → 35.3 (rev 4) → 44.4
(rev 6) — a slow squeeze-limited transient, **not divergence** (increments shrink), settling near
52 W by rev 14. The fix is a **warm-start**: begin each piece at the gas/asperity force balance
(`h_eq = sigma ln(2 P0/|dp|)`), applied automatically when `gas_film_boundary` and
`flat_mixed_lubrication` are both on. Warm-started, the position converges from rev 3 and the
friction from rev 6–8 to ~53 W, so `revolutions=8` is a safe fixed value for the grid study.

## Caveat

The ~53 W bush–vane loss is grid-converged but model-sensitive: it scales roughly linearly with
the asperity scale `FLAT_ASPERITY_PRESSURE_PA` (5–200 MPa), with the assumed 4 MPa bore pressure,
and with the boundary coefficient (0.1). Treat it as order-of-magnitude, not exact.
