# HANDOFF — mochi (rotary/swing compressor physics)

Orientation for continuing the work. The **detailed living roadmap** is
`C:\Users\user\.claude\plans\lazy-puzzling-glade.md`; the **physics of record** is
`PHYSICS.md`. This file is the quick handoff.

---

## 1. 목표 (Goal & current phase)

Incrementally build a physical model of a **swing rotary compressor** (R410A, 30 Hz,
Ø77 bore / Ø68 rotor / e = 4.5 mm). Each rung = **pure-SI module + `PHYSICS.md`
section + `results/` figure + tests**, cross-checked against earlier rungs.

- **Quasi-static phase: COMPLETE** (net DOF = 1, crank angle θ prescribed at constant
  ω, inertia = 0). Chain closes: chamber pressure/volume → gas force/torque →
  indicated work → bearing reaction → drive torque/shaft power + friction + cycle
  performance (leakage, valve).
- **Dynamics (D2) STARTED — the first real time integration is done** (§4.13,
  `rotor_dynamics.py`). The **rotor lateral EOM** `m_r ë = F_gas + F_film(ê,ê̇) +
  m_r ω² O_j` is integrated (implicit BDF) to a steady whirl orbit. Everything under
  it is now **live**: the `mass_properties`/`m_r` inertia, the Ocvirk **ε̇ squeeze
  term**, the numerical Reynolds machinery (`reynolds_1d`, §4.12, matches Ocvirk
  ~1e-4). Coordinate frames all defined (§4.10–4.11): journal `F_j`, bush `F_c/F_f`,
  vane-referenced attitudes `φ_r/φ_b`.
- **Headline dynamic finding:** the **squeeze film limits the eccentricity** — the
  dynamic peak ε ≈ **0.50** (vs the quasi-static **0.71**), min film **7.4 µm** (vs
  4.4). Intrinsic squeeze-lag `τ/T ≈ 0.26` (speed/mass/clearance-independent);
  centrifugal load only ~44 N (~2%). **So §4.9 quasi-static overestimated the peak**;
  the journal runs safer. Reduction check: frozen load → exact quasi-static.
- **Bush multibody (Stage 5) DONE (2026-07):** the two swing-bush pieces are now
  independent bodies in a **9-DOF rotor+bush time integration** (`rotor_bush_dynamics.py`,
  §4.14). Headline: the rotor **gas moment (~22 N·m) is reacted against the fixed vane
  through the bush** (rotor→curved film→piece→flat film→vane), a ~kN load over the 25 mm
  lever that drives the **curved film near contact at peak (min ~0.5 µm, ε_c~0.98)**;
  the rotor bore stays safe (peak ε_j~0.64) and all deviations small (δφ_r~1 mrad).
  **Bush loss of record = dynamic ~0.6 W** (quasi-static 0.2 W kept for comparison only).
  Hydrodynamic-only first cut
  (gas-pressure film bias omitted) — flags the bush as the critical lube site.
- **Rotor–cylinder seal now COUPLED into the 9-DOF by default (2026-08):**
  `seal_contact=True` adds the Hertz line contact + boundary friction to the rotor
  lateral/attitude EOM, so journal + bush + cylinder share load by true stiffnesses.
  **N_c mean ~296N / peak ~1455N (≈2× the free-rotor model), W_r-c ≈ 11 W** — the value
  of record. Attribution: the **bush pulls inward** (−234N); tying the rotor to the fixed
  vane **loads the journal** (+323→+709N outward) which presses the rotor ~1µm deeper →
  N_c doubles. Keep `seal_contact=True` for any **ω-sweep** (seal engagement ∝ ω² via the
  centrifugal). `seal_contact=False` isolates the bush films (reduction checks only).
- **Gas-pressure film BCs (2026-08, `gas_film_boundary=True` opt-in; default keeps the
  validated hydro films):** the 1-D film solvers take chamber/crank gas pressures at the ends
  (`poiseuille_bias_pressure`, superposition) + report throughflow leakage. Piece immersed in
  **bore/recess gas P_bore=4MPa** (user), referenced to bore (divergence theorem). **Both**
  bush films now carry a drop: the **flat vane-film** bore→chamber (IN=suction,
  OUT=compression), and the **curved film** across its **86.5° sealing land** whose ends open
  to the rotor mouth (chamber) and the recess channel (4MPa) — the curved mapping that used
  to be missing (PHYSICS §3.3/4.11). The curved drop gives **+367N/piece away from the vane**,
  largely cancelling the flat film's −400N squeeze.
- **Film cavitation BUG FIXED (2026-08) — supersedes the earlier "×8 bush friction"
  finding:** the clamp used to be switched **off** whenever end pressures were present
  ("gas flooding suppresses cavitation"). At the real state the total **absolute** pressure
  reached **−11.6 MPa** on the diverging half-stroke of the reciprocating pad (−2310N of
  impossible suction). Cavitation constrains the **total** field, so the solvers now solve
  the hydro part unclamped, superpose the gas bias, and clamp the sum at
  `cavitation_pressure_pa` (default 0 = classic Gumbel → **pure-hydro path unchanged**). One
  floor rule for all three films: *no lower than the lowest connected pressure* (degenerates
  to the journal's `p_cav = ambient`). Bush loss in the 3-option regime **91 → 6 W @1800rpm**,
  and the flat film **thickens with speed** instead of thinning.
- **Mixed lubrication = Greenwood-Tripp (2026-08, `mochi.asperity_contact`)** on **both**
  bush films by load sharing, replacing the arbitrary `P0 = 50 MPa` exponential and the hard
  eccentricity clamp. Flat-film friction is split into **boundary vs viscous** on
  `RotorBushOrbit`. **`µ_b` calibration is ruled out** — the coefficient needed to hit any
  target differs **9×** between 1800 and 5400 rpm, so the discrepancy was structural.
  Open: `β, η` need a surface measurement (bush loss spans 66–91 W over the literature range).
- **Thin curved film is a LOAD-PATH consequence, not an artefact:** the attitude EOM closes
  (`M_gas +22.2` vs `M_bush,curved −21.3 N·m`, residual 2% of the reference inertia torque),
  so reacting the gas moment at `l_g`=25mm needs **851N**; an `arc_film` load sweep puts 851N
  at **ε≈0.985 → 0.45µm**, matching the coupled orbit (0.58µm). Adding **Winkler elastic
  compliance** to the arc moved `ε/ε_max` only 1.328→1.330 (the eccentricity is set by
  *geometry* — the piece is pinned to the vane while the groove rides the rotor) and was
  **reverted**. Full record + negative results: `docs/bush_film_revision_2026-08.md`.
- **Still remaining in D2:** 2-D Reynolds; the cavitation floor's exact value (needs an
  oil/refrigerant **solubility** model — bush loss 48.7 W vs 6.2 W across the plausible span);
  GT surface parameters. Rotor EOM, bush multibody, rotor–cylinder contact rung, gas-pressure
  film BCs (both films), and mixed lubrication all exist now.
- **Validation vs Pan et al. 2022 (2026-08) — anchor only, never a fit target.** Different
  machine (valveless swing, 10.34cc) at different pressures, so absolute per-channel losses
  do not transfer; mechanism differences make the efficiency gap unpredictable. What *did*
  transfer: our **journal speed exponent 1.92 vs their 2.00**, and our seal/journal magnitudes
  within ~25% of the size-scaled values — which is what made the 91 W bush channel stand out
  and triggered the review. **No parameter in this repo is fitted to a Pan value.**
- **Bush-film primitives rebuilt (2026-07):** `arc_film` (curved) + `slider_film`
  (flat) now use the **axial-uniform 1-D** reduction (pressure solved along the film
  line × H) via shared `line_reynolds.py`, per the user's spec — NOT the earlier
  axial short-bearing (which was the wrong reduction for the bush L/D>1; kept only for
  the L/D<1 journal). Validated: slider vs Reynolds incline closed form; arc vs slider
  in the shallow limit; manufactured-solution convergence. See PHYSICS §4.11.
- Suite: **416 passed**, `ruff` clean.
- **Rotor-cylinder sealing contact scoped (2026-08):** the swing rotor **slides** on the
  bore (net-zero spin → not a rolling piston), mean |v|≈0.95 m/s (peak 2.0). Boundary
  friction `P=⟨μ N_c |v_slide|⟩` in `rotor_cylinder.py` (μ≈0.1 steel-steel boundary,
  editable; §4.15). Missing input **N_c** (statically indeterminate — same unknown as
  the ~6µm bore-penetration artifact). Quasi-static N_c estimate now in
  `contact_normal_force_n` (journal takes the gas load → N_c~50N mean/190N peak,
  centrifugal 44N floor). **Self-consistent rung DONE (§4.15, `integrate_sealing_contact`):**
  compliant Hertz line contact (Palmgren) in the rotor EOM closes N_c with no free
  param → penetration 6µm→1.4µm (physical), N_c mean ~140N/peak ~780N. Full-EHL check:
  EHL line contact (R_eq=291mm) h_min~0.24µm; at the **Ra 0.3µm design finish** (composite
  RMS σ=0.53µm, Λ=√2·1.25·Ra) Λ~0.45 → **boundary regime**. This standalone
  `integrate_sealing_contact` is a **4-state free-rotor** model → **W_r-c ≈ 6.7 W**; it is
  now a cross-check only. The **coupled 9-DOF value of record is ≈11 W** (see the coupled
  bullet above). Journal-film treatment invalid (c=4.5mm≫thin-film). **So the seal is the
  largest single mechanical loss (≈11W, coupled) — above the journal (8.9W), ≫ bush.**
  `ehl_film_thickness_m` reports Λ.
  Refs: Yanagisawa&Shimizu 1985; mixed-lube review Lubricants 2024; Daikin swing (Purdue/IJR).
- **arc_film validated vs long-bearing analytical (2026-07):** the curved bush film is
  the long-bearing ($L/D>1$) journal limit, so it is cross-checked against the exact
  **infinitely-long (Sommerfeld)** closed form (`long_bearing.py`, §4.11): `arc_film`
  reproduces it to ~1e-7 over ε 0.2–0.9, and carries the attitude invariant
  tan φ = π√(1−ε²)/(2ε). The long-bearing analogue of the §4.12 short-bearing (Ocvirk)
  check — the same analytical-vs-numerical validation at the opposite aspect ratio.

Key numbers (baseline / mouth-aware true geometry): indicated **738 W**, shaft
**747 W** = 738 + bush 0.2 + journal 8.9; peak crank-pin reaction **~2.5 kN**;
η_v **~0.92**; reed-valve overpressure **51 W** (separate performance term, not in the
loads). **Dynamic** journal peak ε **0.50** (quasi-static 0.71). Confirmed masses:
rotor **0.275 kg**, `I_r`=2.11e-4, bush **6.341 g**, shaft 0.379 kg (dormant).

---

## 2. 수정·생성한 파일 (Files)

### `src/mochi/` — physics modules (§ = PHYSICS.md section)
| File | § | Role |
|---|---|---|
| `kinematics.py` | 3.1/4.2 | `RotaryGeometry`, `prescribed_state` (θ→pose), `frequency_hz`/`angular_speed_rad_s`, vane-tip round |
| `rotor_profile.py` | 3.3 | true rotor contour (mouth + lips), `MouthGeometry` |
| `chamber_volume.py` | 3.4 | true-geometry chamber volumes (raster), `SwingBush`, `AxialBands` |
| `chambers.py` | 3.2/3.4 | `port_timed_pressures`, `CycleTrace`/`build_cycle_trace`, seal-over |
| `ports.py` | 3.4 | `characteristic_angles`, `port_open_area_m2` |
| `indicated_work.py` | 3.5 | `W = −∮p dV` → 738 W |
| `thermo_check.py` | 3.5 | CoolProp isentropic triangulation (no measured P-V) |
| `bush_film.py` | 3.6 | swing-bush Couette film friction (~0.2 W) |
| `leakage.py` | 3.7 | orifice leakage + volumetric efficiency |
| `reed_valve.py` | 3.8 | **quasi-static** check valve, overpressure 51 W |
| `gas_force.py` | 4.5 | circular closed-form gas force/torque (fast check, 715 W) |
| `true_gas_force.py` | 4.5 | **mouth-aware** true gas force (contour integral) + torque (virtual work, 738 W); `peak_rotor_force_n` (~2.5 kN) — the **headline** |
| `bearing_load.py` | 4.6 | `R_j = −F_gas`, drive torque, `shaft_work_j` → **747 W** (baseline) |
| `journal_bearing.py` | 4.7 | Petroff journal friction 8.9 W (CAD `r_j = 14.2 mm`) |
| `mass_properties.py` | 4.8 | **D2 groundwork** — density-editable rotor/bush mass + polar inertia (**data only; unused by any solver yet**) |
| `ocvirk_bearing.py` | 4.9 | **D2 groundwork (quasi-static)** — L1 closed-form short-bearing (Ocvirk) journal force `F(ε, ε̇, Ω)`, **not a PDE**; `eccentricity_cycle` = **steady-ε force balance** (ε̇=0), running `ε(θ)` from `R_j(θ)`, eccentric friction, Petroff reduction |
| `reynolds_1d.py` | 4.10–4.12 | **D2** — 1-D finite-width (short-bearing) **numerical Reynolds** solver (`solve_short_bearing_1d`): tridiagonal axial FD + Gumbel cavitation, **validates vs Ocvirk to ~1e-4**. First PDE-style film solver; machinery for the later 2-D |
| `long_bearing.py` | 4.11 | **D2 — infinitely-long (Sommerfeld) journal analytical** (`long_bearing_load`, `sommerfeld_pressure`): the closed-form long-bearing ($L/D>1$) benchmark for the curved bush film. `arc_film` matches it to ~1e-7 (ε 0.2–0.9); attitude invariant tan φ=π√(1−ε²)/(2ε). The long-bearing twin of the Ocvirk short-bearing check |
| `rotor_dynamics.py` | 4.13 | **NEW (D2) — the rotor lateral EOM + time integration** (`integrate_rotor_orbit` → `RotorOrbit`): `m_r ë = F_gas + F_film(ê,ê̇) + m_r ω² O_j`, implicit BDF (`scipy.solve_ivp`) → steady whirl orbit. **First real mechanical time integration.** Dynamic peak ε 0.50 < quasi-static 0.71 (squeeze lag). `ROTOR_MASS_KG`=0.275 editable |
| `rotor_bush_dynamics.py` | 4.14 | **NEW (D2, Stage 5) — 9-DOF rotor + two-piece swing-bush multibody** (`integrate_rotor_bush_orbit` → `RotorBushOrbit`): 18-state stiff BDF; both bush pieces are bodies coupled by the curved (`arc_film`) + flat (`slider_film`) films. Gas moment reacted through the bush → curved film near contact at peak (~0.5 µm). Confirmed `I_r`=2.11e-4, `m_p`=6.341 g. Stiff — loose velocity `atol` lets BDF step over the squeeze mode. **`seal_contact=True` (default)** couples the rotor–cylinder Hertz contact into the EOM (N_c ~296N, W_r-c ~11W; keep on for ω-sweeps) |
| `asperity_contact.py` | 4.11 | **NEW (2026-08)** — Greenwood-Tripp elastic asperity contact `p_asp(h)` from measurable surface parameters (σ, β, η, E'), replacing an arbitrary pressure scale. Load-shares with **both** bush films; `AsperityParams` carries the literature ground-steel defaults (ηβσ≈0.05) pending a surface measurement |
| `bush_outline.py` | 3.3 | **NEW (2026-08)** — the single source for the swing-bush piece and vane outlines in mm, shared by the results figures and `bush_gui` (was copy-pasted in three places). `close_bore_arc` selects the vane's top bore arc |
| `bush_gui.py` | 4.14 | **NEW (2026-08)** — coupled-orbit viewer: true-scale rotor-mouth panel + ×60 clearance inset, animated over one revolution. GIF export or Tk window with a crank-angle slider (`mochi bush-gui`). Orbit is cached to `.npz` since it costs minutes |
| `tecplot.py` | — | **NEW (2026-08)** — Tecplot ASCII writer (1-D point zones, 2-D ordered zones, and **transient strands** with STRANDID/SOLUTIONTIME). Backs the `results/data/` raw export from `generate_results.py --data`; see `docs/data_export.md` |
| `gui.py`, `cli.py`, `__main__.py` | — | visualization / plumbing (`cli.py` also fronts `bush-gui`) |

### tests/ — `test_*.py` per module (incl. `test_reynolds_1d.py` 15, new `test_rotor_dynamics.py` 9)
### `scripts/generate_results.py` — every figure under `results/`; new coordinate figs `geometry/{journal_film_coordinates, journal_film_axial, bush_film_coordinates, bush_attitude_reference}.png` + `bearing_load/{reynolds_1d_validation, rotor_orbit}.png`
### `PHYSICS.md` — §3.5–3.8, §4.5–**4.13**, §9 open items
### `docs/HANDOFF.md` — this file

---

## 3. 내린 결정 (Key decisions)

- **True geometry is the headline** for gas force/torque/bearing (mouth cavity lowers
  the rotor load ~20% vs the circular disc; torque via **virtual work** = 738 W
  exactly). Circular closed form kept as the fast analytic cross-check.
- **Baseline (reed-valve-free) unification**: all mechanical **loads** use the ideal
  indicated 738 W; the reed-valve overpressure (51 W → valve-aware 789) is a
  **separate §3.8 performance term, NOT propagated into the loads** (shown only as a
  labelled overlay on the P-V diagram). Shaft power = **747 W**.
- **Rotor–cylinder contact friction: omitted** — it is a rolling/sliding line contact
  at eccentricity ε ≈ 1 (R_c−R_r = 4.5 mm = e), so a concentric Petroff is invalid;
  it needs the indeterminate contact force (EHL/compliance) → deferred to D2.
- **CAD journal** `r_j = 14.2 mm` (was assumed 16); `c_j = 15 µm` still assumed.
- **φ notation** (to avoid clashing with the suction-seal angle): **rotor swing =
  `\phi` (straight ϕ)**, **suction seal angle = `\varphi` (curly φ)** — in figures and
  `PHYSICS.md`.
- **MBD scope decisions** (for the next phase):
  - **ω is PRESCRIBED CONSTANT** (operating condition is a drive frequency in Hz) →
    **D3 (free θ / speed ripple) is DEFERRED**; no motor torque / crank `I_eff` needed.
  - **D1 (reed valve) kept quasi-static**; the dynamics **start at D2**.
  - **D2 = "rotordynamics at fixed drive speed"** (rotor small motion in bearing
    clearances under the quasi-static gas load). Target fidelity **L2 (full Reynolds)**
    for the films (L1 Ocvirk as a build-up/fallback); rotor–cylinder = separate
    Hertz/EHL contact.
  - Mass/inertia + reed data sourced **hybrid** (CAD/measured where available, else
    assumed). **Density is EDITABLE** (rotor/bush material not confirmed).

---

## 4. 남은 작업 (Remaining work)

### D2 — STARTED (the rotor EOM is integrated; multi-body picture still partial)
The dynamics rung proper — a rotor equation of motion integrated in time — is **now
done** (`rotor_dynamics.py`, §4.13, see below). Its quasi-static building blocks:
- **✓ Groundwork — `mass_properties.py`**: density-editable rotor/bush mass + inertia.
  **Data only** — no module imports it; it feeds no equation of motion yet.
- **✓ Groundwork — `ocvirk_bearing.py`** (§4.9): **L1 closed-form short-bearing
  (Ocvirk) journal force** `F(ε, ε̇, Ω)` — π-film, **not a PDE** (L2 = full Reynolds
  is deferred). Radial/tangential components, textbook load capacity + attitude
  angle (machine-precision tests), `equilibrium_eccentricity_ratio` (inverts load→ε),
  eccentric friction `T_petroff/√(1−ε²)` (→ Petroff as ε→0, the reduction check).
  `eccentricity_cycle` is a **steady-eccentricity force balance per crank angle
  (ε̇ = 0), NOT a time integration**: it balances the film against the true `R_j(θ)`
  → **running ε 0.13–0.71 (peak 0.71 @ 226°), min film ~4.4 µm, journal loss
  9.2→10.8 W (~18%)** — quantifying the §4.7 concentric-Petroff underestimate. The
  ε̇ squeeze term exists in `short_bearing_force` but is **never driven nonzero**
  outside tests.
- **✓ Groundwork — coordinate frames + 1-D Reynolds solver** (§4.10–4.12):
  - **Frames defined** (figures only, no dynamic code): crank-pin journal `F_j`
    (Cartesian `ê=O_b−O_j`, symmetric 4 MPa BC, offset-invariant), swing-bush curved
    `F_c` (partial journal) + flat `F_f` (slider), **vane-referenced attitudes**
    `φ_r=φ−90°` (rotor, prescribed) / `φ_b` (bush, now a small DOF; `φ_b=0` recovers
    §3.6/§4.7). Three-centre convention nailed: O (shaft axis) / O_j (crank-pin) /
    O_b=O_r (rotor bore); crank-throw `e`=4.5 mm ≠ bearing ecc `|ê|`≤15 µm.
  - **`reynolds_1d.py`** (§4.12): 1-D finite-width **numerical Reynolds** solve
    (tridiagonal axial FD + Gumbel cavitation) → **matches Ocvirk to ~1e-4** (static);
    squeeze differs ~2% (π-film superposition vs true cavitation region — honest).
    The **numerical machinery** the 2-D solver reuses. Still not the coupled EOM.
- **✓ DONE — Rotor lateral EOM + time integration** (§4.13, `rotor_dynamics.py`):
  `m_r ë = F_gas + F_film(ê,ê̇) + m_r ω² O_j` integrated with implicit BDF
  (`scipy.solve_ivp`) to a steady whirl orbit. **First real mechanical time
  integration.** The ε̇ squeeze term + m_r inertia are now live. **Finding:** the
  squeeze film **attenuates+lags** the swing → dynamic peak ε **0.50 < quasi-static
  0.71**, min film **7.4 µm**, whirl orbit a small bounded loop; centrifugal ~44 N
  (~2%). Validated: static balance exact, frozen-load → quasi-static, periodic.
  Figure `bearing_load/rotor_orbit.png`; tests (9). Uses `c_j`=15 µm, `m_r`=0.275.
- **NEXT — remaining D2, in order:**
  1. **c_j sweep** {7.5, 15, 30} µm on the dynamic orbit (the peak-ε sensitivity;
     data is user-confirmable, currently 15 µm).
  2. **Rotor–cylinder contact** (Hertz/EHL) → determinate contact normal force `N_c`
     (resolves the §4.6 indeterminacy) + contact friction. Needs contact stiffness+μ
     (not in the model — see §6).
  3. **✓ DONE — Bush as a moving body** (Stage 5, §4.14, `rotor_bush_dynamics.py`):
     the two pieces are independent bodies (curved `F_c` + flat `F_f` films, attitudes
     `φ_k`), 9-DOF stiff BDF. Next on this thread: add the **gas-pressure film bias**
     (§9.3) and an EHL/contact treatment for the peak-load near-contact curved film.
  4. **Upgrade films to L2 (2-D Reynolds)** — reuse the `reynolds_1d` machinery on an
     (α,ζ) grid with tilt; the biggest numerical piece.
  - **Data still needed** (see the new §6 inventory): rotor/bush **densities** and
    bearing **clearances** are user-confirmable (currently assumed 7200/7850,
    15/10/8.5 µm); rotor–cylinder **contact stiffness + μ** and reed **m/k/c** do
    **not exist in the model at all** and must be sourced before those rungs.
  - **Validation target**: quasi-static → dynamic **reduction check** (as inertia→0 /
    stiffness→∞ / Reynolds→concentric, must recover the quasi-static results — the
    §4.9 `eccentricity_cycle` result is the target the dynamic film must reproduce).

### Deferred / parked
- **D1** transient reed valve (spring-mass-damper lift ODE; back-flow/flutter) — only
  when the valve overpressure must feed the loads.
- **D3** free θ (speed ripple) — only if cycle speed ripple matters (needs `I_eff` +
  drive torque).
- **D4** dynamic gas state (valve↔pressure↔motion coupling; heat transfer → real
  indicated +15–40%).
- Non-dynamics: bush-film temperature-dependent viscosity; short-bearing already
  folds into D2.

### PHYSICS.md §9 open items still unchecked
- Confirm the prototype frame against CAD; recess/bush pressure BCs (fixed 4.0 MPa);
  oil-film leakage across recess/bush seals; reed dynamics/back-flow + true port
  size/Cd vs measurement.

---

## 5. 검증 명령어 (Verification)

```bash
# Full gate (from repo root)
ruff format --check . && ruff check . && pytest        # → 335 passed, ruff clean
#  (mypy still fails on a pre-existing NumPy-stub / py3.11 issue — unrelated.)

# Regenerate every figure (no glyph warnings expected)
python scripts/generate_results.py

# Spot-checks (headline numbers)
python -c "from mochi.bearing_load import shaft_work_j; from mochi.kinematics import RotaryGeometry; w=shaft_work_j(RotaryGeometry.default()); print(w.shaft_power_w, w.journal_friction_power_w, w.peak_journal_force_n)"
#   → shaft ~747 W, journal ~8.9 W, peak ~2.5 kN

python -c "from mochi.true_gas_force import true_gas_torque_work_j, peak_rotor_force_n; from mochi.kinematics import RotaryGeometry; g=RotaryGeometry.default(); w=true_gas_torque_work_j(g); print(w.power_w, w.circular_power_w, peak_rotor_force_n(g))"
#   → true torque ~738 W (= indicated), circular ~715 W, peak ~2.5 kN

python -c "from mochi.mass_properties import rotor_mass_properties, bush_mass_properties; from mochi.kinematics import RotaryGeometry; g=RotaryGeometry.default(); r=rotor_mass_properties(g); b=bush_mass_properties(g); print(round(r.mass_kg,4), r.inertia_kg_m2, round(b.mass_kg,4), b.inertia_kg_m2)"
#   → rotor ~0.42 kg / ~2.8e-4 kg·m² ; bush ~0.013 kg / ~5.5e-7 kg·m²  (edit ROTOR_/BUSH_DENSITY_KG_M3 to rescale)

python -c "from mochi.leakage import leaky_cycle; from mochi.kinematics import RotaryGeometry; c=leaky_cycle(RotaryGeometry.default()); print(round(c.capped_peak_pa/1e6,1), round(c.volumetric_efficiency,3))"
#   → ~9.7 MPa, η_v ~0.92

python -c "from mochi.ocvirk_bearing import eccentricity_cycle; from mochi.kinematics import RotaryGeometry; c=eccentricity_cycle(RotaryGeometry.default()); print(round(c.peak_eccentricity_ratio,3), round(c.minimum_film_thickness_m*1e6,2), round(c.mean_friction_power_w,1), round(c.petroff_mean_friction_power_w,1))"
#   → peak ε ~0.71, min film ~4.36 µm, eccentric 10.8 W vs Petroff 9.2 W  (§4.9; slow — cold trace ~75 s)

python -c "from mochi.reynolds_1d import solve_short_bearing_1d; from mochi.ocvirk_bearing import short_bearing_force; n=solve_short_bearing_1d(0.7,95.0); o=short_bearing_force(0.7,95.0); print(round(n.magnitude_n,2), round(o.magnitude_n,2), f'{abs(n.magnitude_n-o.magnitude_n)/o.magnitude_n:.1e}')"
#   → 1D Reynolds ≈ Ocvirk (§4.12), rel err ~1e-4

python -c "from mochi.rotor_dynamics import integrate_rotor_orbit; from mochi.kinematics import RotaryGeometry; o=integrate_rotor_orbit(RotaryGeometry.default()); print(round(o.peak_eccentricity_ratio,3), round(o.quasi_static_peak_eccentricity_ratio,3), round(o.minimum_film_thickness_m*1e6,2), round(o.centrifugal_load_n,1))"
#   → dynamic peak ε ~0.50 < quasi-static ~0.71, min film ~7.4 µm, centrifugal ~44 N  (§4.13; slow — cold trace ~75 s)
```

---

## 6. 물성치 확보 상태 (Material properties — assumed vs confirmed)

Every physical input the model **assumes** vs what is **confirmed** (CAD / supplied /
CoolProp). The dynamics rungs' fidelity is bounded by group A/B; **fill A with real
values, and B must be sourced before the contact / reed / L2 rungs.**

### A. Confirmable from the real parts — currently ASSUMED, replace when available
| Property | Assumed value | Where (constant) |
|---|---|---|
| Rotor material / density ρ_r | 7200 kg/m³ (grey cast iron) | `mass_properties.py` `ROTOR_DENSITY_KG_M3` |
| Bush material / density ρ_b | 7850 kg/m³ (steel) | `mass_properties.py` `BUSH_DENSITY_KG_M3` |
| Crank-pin journal radial clearance c_j | 15 µm (not drawn — nominal fit) | `journal_bearing.py` `JOURNAL_CLEARANCE_M` |
| Swing-bush film clearance | 10 µm (from the piece-shift geometry) | `bush_film.py` (SwingBush `piece_shift`) |
| Axial running clearance | 8.5 µm / side | `bush_film.py` `AXIAL_CLEARANCE_M` |
| Oil viscosity μ (operating) | 0.010 Pa·s (POE VG68 @ ~80 °C, dilution neglected) | `bush_film.py` `LUBRICANT_VISCOSITY_PA_S` |

Editing the density constant (or passing `density_kg_m3`) rescales mass + inertia
linearly; `c_j` sets the running film margin of §4.9 directly.

### B. Not in the model at all — must be sourced before the dependent dynamics rung
| Property | Needed by | Status |
|---|---|---|
| Reed-valve equivalent mass / stiffness / damping / max lift | transient reed valve (D1) | **absent** — §3.8 is a quasi-static check valve (only preload 5 %, Cd 0.6, geometric port area) |
| Rotor–cylinder contact stiffness (Hertz/EHL) + friction coeff. μ | contact rung (resolves §4.6 indeterminacy) | **absent** — contact force omitted, no coefficient |
| Oil cavitation pressure | L2 (full Reynolds) film | **absent** — only needed once the PDE film is built |

### C. Confirmed — CAD / supplied / CoolProp (do not change)
Journal `r_j = 14.2` mm, `L_j = 21` mm (CAD); supplied ports 0.82 / 3.24 MPa;
polytropic `n = 1.07`, R410A `ρ_suc = 31.4`, nozzle `k = 1.10` (CoolProp HEOS);
geometry Ø77 / Ø68 / e = 4.5 mm + port timing + bush/mouth dims (supplied 2026-07).

> **Note on scope:** §4.8/§4.9 are **quasi-static D2 groundwork** (a force law + inertia
> data), **not dynamics** — no equation of motion or time integration exists on the
> mechanical side yet (see §4).
