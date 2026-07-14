# Physics and numerical methods

This document is the contract between the physical model and its
implementation. It starts deliberately conservative: unconfirmed compressor
geometry and force laws are listed as open decisions rather than presented as
facts.

## 1. Problem statement and initial scope

The project will become a broader rotary-compressor simulation program. The
current accepted implementation is only a supporting prescribed-motion
visualization for checking mechanism geometry and animation conventions.

Currently in scope:

- a fixed cylinder and a changing-length vane attached at the top center;
- a prescribed eccentric orbit for the rotor center;
- a rotor-fixed circular cutout constrained to the vane centerline;
- a white clearance slot slightly wider than the vane;
- physical crank speed separated from slowed display speed.

Not implemented by this visualization:

- pressure, flow, force, contact, friction, torque, or deformation;
- initial-value dynamics or numerical time integration;
- bearing, lubrication, thermal, or fluid-structure models;
- the main rotary-compressor solver and its production outputs.

These belong to the broader simulation only after their equations and
validation cases are reviewed.

## 2. Coordinates, signs, and units

All internal quantities use SI units: metre, kilogram, second, newton, pascal,
joule, kelvin, and radian. Input/output adapters may display other units only
through explicit conversion.

The visualization frame is fixed to the cylinder and defined in Section 3.
No unconstrained generalized body state or force convention is implemented.
Any later dynamics module must define its own generalized coordinates, force-on-
body convention, and mapping to this visualization frame before integration.

Angles are in radians internally. Each angular input must state its zero
direction and positive sense. Pressure models must state whether pressure is
absolute or gauge and which reference pressure is used.

## 3. Prescribed rotary mechanism

The test GUI uses a fixed cylinder-centered frame. The origin is the cylinder
center, positive $y$ points toward the top vane, and positive $x$ points toward
the inlet side. The crank angle $\theta=0$ places the rotor at the top of its
orbit, $(x_r,y_r)=(0,e)$; increasing $\theta$ is clockwise.

The supplied geometry is:

| Symbol | Meaning | Value |
|---|---|---:|
| $D_c$ | cylinder inside diameter | 77.0 mm |
| $D_r$ | main rotor outside diameter | 68.0 mm |
| $e$ | rotor-center eccentricity | 4.5 mm |
| $r_h$ | circular-cutout radius | 8.0 mm |
| $L$ | rotor-center to cutout-center distance | 25.0 mm |
| $w_v$ | vane width | 8.0 mm |
| $a_v$ | rotor-center to vane-tip distance at the top position | 9.0 mm |
| $f$ | physical rotation frequency | 30 Hz = 1800 rpm |

The eccentricity is the radial difference,

$$
e = \frac{D_c-D_r}{2} = 4.5\ {\rm mm},
$$

so the main rotor remains internally tangent to the cylinder in this ideal
geometry. Its prescribed center is

$$
x_r=e\sin\theta,\qquad y_r=e\cos\theta.
$$

The circular-cutout center $\mathbf{h}=(0,y_h)$ is fixed in the rotor a distance
$L$ from the rotor center and constrained to the stationary vane centerline. Selecting
the upper assembly branch gives

$$
y_h=y_r+\sqrt{L^2-x_r^2}.
$$

The rotor orientation is the direction from the rotor center to the cutout center,

$$
\phi=\operatorname{atan2}(y_h-y_r,-x_r).
$$

The cutout fits inside the 34 mm rotor radius because $L+r_h=33$ mm. A radial
slot opens the circular cutout through the rotor OD and stops at the circular
boundary; no rectangular or square cutout continues toward the rotor center.
Its default displayed width is $w_s=w_v+1$ mm. At the two rotor lips, tangent
fillets with nominal radius $r_f=1.5$ mm round the transition from the rotor OD
to the slot. The displayed fillet radius is reduced automatically only if
edited inputs leave insufficient material.

Let $\mathbf{r}=(x_r,y_r)$ and define the rotor-fixed axial and transverse unit
vectors

$$
\mathbf{u}=\frac{\mathbf{h}-\mathbf{r}}{L},\qquad
\mathbf{n}=(-u_y,u_x).
$$

For the prescribed vane animation only, an invisible rotor-fixed reference line
passes through $\mathbf{c}=\mathbf{r}+a_v\mathbf{u}$ in the
$\mathbf{n}$ direction. This reference is not a cutout boundary. The vertical
vane remains on $x=0$, so its lower tip is the intersection

$$
\lambda=-\frac{c_x}{n_x},\qquad
\mathbf{p}_{tip}=\mathbf{c}+\lambda\mathbf{n}.
$$

At $\theta=0$, $\lambda=0$ and the tip is exactly $a_v=9$ mm from the
rotor center. At other angles the intersection moves, and the displayed
centerline vane length is

$$
\ell_v=\frac{D_c}{2}-y_{tip}.
$$

The light-gray rotor is one parametric polygon built from the rotor OD, two
fillet arcs, the slot sides, and the remaining circular-cutout circumference.
The majority cylinder arc and the right, bottom, and left vane sides are one
continuous line with no fill. Its interior is transparent, so the vane can
overlap the rotor below the circular cutout; the surrounding canvas remains
white. No white shapes are overpainted to manufacture a cutout.

At constant physical frequency,

$$
\theta(t)=\theta_0+2\pi f t.
$$

The display slow factor $s$ changes wall-clock animation only:
$\theta_{display}(t_{wall})=\theta_0+2\pi f s t_{wall}$. It must never alter
physical time, velocity, force, or result metadata.

The cylinder center `C`, rotor center `R`, and dashed circular locus of radius
$e$ are display guides only; they do not add geometry or a physical boundary.

Port markers are display references rather than flow boundaries. Angles are
measured from positive $y$ toward positive $x$: the inlet is at $+30^\circ$
(top right) and the approximate symmetric outlet is at $-30^\circ$ (top left).

This prototype is prescribed kinematics. The cutout centerline constraint,
ideal tangency, and displayed slot clearance are not contact-force or
lubrication models.

## 4. Future dynamics placeholder (not implemented)

The visualization evaluates no equation of motion. A future dynamics module
may introduce generalized position $\mathbf{q}$, velocity
$\mathbf{v}=\dot{\mathbf{q}}$, and mass matrix $\mathbf{M}$:

$$
\mathbf{M}\ddot{\mathbf{q}} = \mathbf{F}_{total}.
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

If dynamic rotation is added, define an orientation $\phi$, moment of inertia
$I$, and a moment balance $I\ddot{\phi} = M_{total}$ with the same explicit
sign discipline. Do not infer rotation from planar force data alone.

## 5. Kinematics and force interfaces

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

## 6. Initial-value formulation

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

## 7. Required validation ladder

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

## 8. Case and result contract

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

## 9. Open model decisions

The following questions block a physical production solver:

- [ ] Define the future dynamic bodies and their interfaces.
- [x] Define the prototype global frame, angle origin, and rotation direction.
- [ ] Confirm that prototype frame against the compressor CAD convention.
- [ ] Decide the required degrees of freedom and constraints.
- [x] Document the supplied cylinder, rotor, circular-cutout, and vane dimensions.
- [ ] Confirm the invisible vane-tip reference relation and remaining clearances against CAD.
- [ ] Identify all applied, gas, contact, friction, and fluid forces.
- [ ] Choose whether time or crank angle is the primary independent variable.
- [ ] Select the contact/friction treatment, if needed.
- [ ] Establish at least one trusted compressor reference case.
- [ ] Set accuracy and conservation acceptance criteria.

Record resolved decisions here with equations and evidence; do not merely check
the box.

## 10. Model-change procedure

Any pull request that changes physics or numerics must update this document in
the same change. State the old and new equation, assumptions, expected effect,
and validation evidence. Update `CHANGES.md`, and add or revise tests that would
fail under the wrong sign, unit, or model.
