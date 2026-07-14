# Physics and numerical methods

This document is the contract between the physical model and its
implementation. It starts deliberately conservative: unconfirmed compressor
geometry and force laws are listed as open decisions rather than presented as
facts.

## 1. Problem statement and initial scope

The project will predict the motion of a bushing in a scroll compressor from a
specified initial state and a set of explicit force contributions. The first
accepted model is intended to cover planar rigid-body translation. Additional
degrees of freedom and couplings enter only after their equations and
validation cases are reviewed.

Initially in scope:

- bushing-center position and velocity in one defined planar frame;
- time- or crank-angle-dependent kinematics supplied by the case definition;
- modular force contributions acting on the bushing;
- an initial-value time integration with visible tolerances and failure state;
- histories of state, forces, and solver diagnostics.

Not yet accepted by default:

- bushing rotation or axial motion;
- contact, impact, clearance, or friction laws;
- hydrodynamic or mixed-lubrication film forces;
- gas-pressure, orbiting-scroll, bearing, or drive-force closures;
- thermal deformation, wear, or fluid-structure coupling.

These may be important, but their absence here means ?to be defined,? not ?to
be neglected in the real compressor.?

## 2. Coordinates, signs, and units

All internal quantities use SI units: metre, kilogram, second, newton, pascal,
joule, kelvin, and radian. Input/output adapters may display other units only
through explicit conversion.

The global planar frame must be tied to compressor geometry before the first
physical force model is merged. Until then, use the following mathematical
symbols without assigning a machine-specific direction:

$$
\mathbf{q}(t) =
\begin{bmatrix}x_b(t) \\ y_b(t)\end{bmatrix},\qquad
\mathbf{v}(t) = \dot{\mathbf{q}}(t).
$$

Here, $\mathbf{q}$ is the bushing-center displacement in the chosen inertial
frame. A positive force component accelerates the bushing in the corresponding
positive coordinate direction. Every force name and equation must specify
?force on the bushing?; reaction forces on other components have the opposite
sign.

Angles are in radians internally. Each angular input must state its zero
direction and positive sense. Pressure models must state whether pressure is
absolute or gauge and which reference pressure is used.

## 3. Baseline equation of motion

For planar translation with constant bushing mass $m_b$:

$$
m_b\ddot{\mathbf{q}} = \mathbf{F}_{total}.
$$

The total is assembled from named contributions:

$$
\mathbf{F}_{total} =
\mathbf{F}_{drive} +
\mathbf{F}_{gas} +
\mathbf{F}_{contact} +
\mathbf{F}_{fluid} +
\mathbf{F}_{support} +
\mathbf{F}_{other}.
$$

Only accepted, configured terms are evaluated. An absent model contributes
zero and must be reported in case metadata; it must not be approximated
silently.

If a linear support is later justified, a possible constitutive term is

$$
\mathbf{F}_{support} =
-\mathbf{C}\mathbf{v}
-\mathbf{K}(\mathbf{q}-\mathbf{q}_{ref}),
$$

where $\mathbf{C}$ has units N s/m and $\mathbf{K}$ has units N/m. This is an
interface example, not yet an accepted model of the compressor.

If bushing rotation is added, define an orientation $\phi_b$, moment of inertia
$I_b$, and a moment balance $I_b\ddot{\phi}_b = M_{total}$ with the same explicit
sign discipline. Do not infer rotation from planar force data alone.

## 4. Kinematics and force interfaces

A compressor case must provide the geometry and prescribed driver state needed
by its force laws. If crank speed is prescribed,

$$
\theta_c(t) = \theta_{c,0} + \int_0^t \omega_c(\tau)\,d\tau.
$$

Constant speed is a special case, not a universal assumption. The crank-angle
origin and positive rotation must be documented with the first case.

Each force contribution should be a deterministic function with an interface
equivalent to

$$
\mathbf{F}_i =
\mathbf{F}_i(t,\mathbf{q},\mathbf{v},\text{case data}).
$$

It must document:

- the component exerting the force and the body receiving it;
- the frame in which inputs and outputs are expressed;
- required geometry and parameter units;
- validity range and singular or discontinuous conditions;
- source, derivation, or reference data;
- an independent test.

For unilateral contact, the eventual model must at least respect the gap and
normal-reaction conditions

$$
g(\mathbf{q}) \ge 0,\qquad \lambda_n \ge 0,\qquad
\lambda_n g(\mathbf{q}) = 0,
$$

before a penalty, complementarity, event-based, or other numerical treatment
is chosen. Friction requires a separate tangential law and sign convention.

## 5. Initial-value formulation

Define the first-order state

$$
\mathbf{y} =
\begin{bmatrix}\mathbf{q} \\ \mathbf{v}\end{bmatrix},\qquad
\dot{\mathbf{y}} =
\begin{bmatrix}
\mathbf{v} \\
\mathbf{F}_{total}(t,\mathbf{q},\mathbf{v})/m_b
\end{bmatrix}.
$$

For smooth forces, the first reference implementation may use SciPy's adaptive
`solve_ivp` interface with an explicit method and reported relative/absolute
tolerances. Output sampling must be independent of internal solver steps.

This choice must be reconsidered when discontinuous contact, stiffness, or
algebraic constraints are introduced. An adaptive solver reporting success is
not by itself evidence that the physical solution is correct.

State variables with different magnitudes need component-appropriate absolute
tolerances or documented nondimensionalization. A case must record the method,
tolerances, step limits, number of accepted/rejected steps when available, and
termination reason.

## 6. Required validation ladder

Implementations should progress through these tests in order:

1. **Zero force:** velocity is constant and position is linear in time.
2. **Constant force:** acceleration is constant and matches the analytic
   quadratic trajectory.
3. **Linear oscillator:** if a support term is implemented, compare frequency
   and damping decay with the analytic solution.
4. **Frame/sign cases:** simple axis-aligned loads produce motion in the stated
   positive direction; action/reaction pairs have opposite signs.
5. **Step/tolerance refinement:** quantities of interest approach a stable
   value as error tolerances tighten or maximum step decreases.
6. **Energy/work balance:** for conservative cases, total energy remains within
   the stated numerical error; for dissipative cases, the loss has the correct
   sign and magnitude.
7. **Compressor reference case:** compare against independently calculated,
   experimental, or trusted published results with provenance and acceptance
   limits.

Regression tolerances must reflect a measured numerical error. They must not be
chosen solely to accept the current output.

## 7. Case and result contract

Every reproducible case should eventually record:

- model version or Git commit;
- geometry and material parameters with units;
- active and inactive force models;
- initial state and simulated interval;
- prescribed speed/load histories and their interpolation rule;
- numerical method, tolerances, and step limits;
- requested outputs and output times.

Every result should include time, position, velocity, total force, each active
force contribution, and solver status. Results from failed or incomplete solves
must remain clearly marked and must not look like successful final states.

## 8. Open model decisions

The following questions block a physical production solver:

- [ ] Define the physical bushing and its interfaces to neighboring components.
- [ ] Fix the global coordinate frame, angle origin, and rotation direction.
- [ ] Decide the required degrees of freedom and constraints.
- [ ] Document compressor geometry and nominal clearances.
- [ ] Identify all applied, gas, contact, friction, and fluid forces.
- [ ] Choose whether time or crank angle is the primary independent variable.
- [ ] Select the contact/friction treatment, if needed.
- [ ] Establish at least one trusted compressor reference case.
- [ ] Set accuracy and conservation acceptance criteria.

Record resolved decisions here with equations and evidence; do not merely check
the box.

## 9. Model-change procedure

Any pull request that changes physics or numerics must update this document in
the same change. State the old and new equation, assumptions, expected effect,
and validation evidence. Update `CHANGES.md`, and add or revise tests that would
fail under the wrong sign, unit, or model.
