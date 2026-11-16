# RVE Studio 1.0 — Theory and Modelling Guidance

This document states every equation RVE Studio solves, the assumptions behind each wall model, the limits within which
the one-dimensional wall model has been verified against shell finite elements, and how to decide whether a finite
element (FE) result is mesh-converged. It is the reference for choosing a wall model and for judging a result.
The step-by-step use of the program is in the [User Guide](USER_GUIDE.md).

**Source model.** J. Sriharan, M. Dias, D. Fernando, S. Adhikari (2024), *A unified approach for the prediction of
the effective properties of laminated composite cellular core*, Engineering Structures 299, 117106,
doi:[10.1016/j.engstruct.2023.117106](https://doi.org/10.1016/j.engstruct.2023.117106). RVE Studio reproduces the 84
reference values of that paper within 0.3 % and extends it with shear-flexible walls, a calibrated
twist-corrected wall, a reliability rating, periodic shell FE validation with mesh-convergence verdicts, design
databases and wall-buckling screening.

**How to cite.** RVE Studio builds on earlier papers; please cite the RVE Studio paper together with them, and the software
release through its Zenodo DOI (`CITATION.cff`):

1. J. Sriharan et al. RVE Studio: open-source strain-energy homogenisation of laminated-composite cellular cores —
   verification of the full stiffness tensor, a bend–twist reliability limit and design tools. *In preparation.*
2. J. Sriharan, M. Dias, S. Adhikari, D. Fernando (2025). A semi-analytical model for predicting the shear buckling of
   laminated composite honeycomb cores in sandwich panels. *Composite Structures* 351, 118629.
   doi:[10.1016/j.compstruct.2024.118629](https://doi.org/10.1016/j.compstruct.2024.118629)
3. J. Sriharan, M. Dias, D. Fernando, S. Adhikari (2024). A unified approach for the prediction of the effective
   properties of laminated composite cellular core. *Engineering Structures* 299, 117106.
   doi:[10.1016/j.engstruct.2023.117106](https://doi.org/10.1016/j.engstruct.2023.117106)

---

## Contents

1. [Notation and conventions](#1-notation-and-conventions)
2. [Ply and laminate](#2-ply-and-laminate)
3. [Equivalent transverse shear stiffness of a wall](#3-equivalent-transverse-shear-stiffness-of-a-wall)
4. [The strain-energy wall element](#4-the-strain-energy-wall-element)
5. [Wall bending options](#5-wall-bending-options)
6. [Periodic RVE and homogenisation](#6-periodic-rve-and-homogenisation)
7. [Effective properties](#7-effective-properties)
8. [Physical properties](#8-physical-properties)
9. [Reliability of the wall model](#9-reliability-of-the-wall-model)
10. [Choosing a wall model — guidance](#10-choosing-a-wall-model--guidance)
11. [Benchmark designs](#11-benchmark-designs)
12. [Shell FE validation](#12-shell-fe-validation)
13. [Mesh convergence: is my FE result converged?](#13-mesh-convergence-is-my-fe-result-converged)
14. [Wall load recovery](#14-wall-load-recovery)
15. [Local wall buckling](#15-local-wall-buckling)
16. [Design databases and post-processing](#16-design-databases-and-post-processing)
17. [Limitations](#17-limitations)
18. [References](#18-references)

---

## 1. Notation and conventions

| Symbol | Meaning |
|---|---|
| $x_1, x_2$ | in-plane axes of the core (cell plane); $x_3$ the core height direction |
| $\mathbf a_1, \mathbf a_2$ | periodic lattice vectors of the RVE; RVE area $A_{RVE} = \lvert \mathbf a_1 \times \mathbf a_2 \rvert$ |
| $h$ | core height; RVE volume $V = A_{RVE}\,h$ |
| $\bar 1, \bar 2, \bar 3$ | wall-local axes: $\bar 1$ along the wall from Node I to Node J, $\bar 3 = x_3$, $\bar 2 = \bar 3 \times \bar 1$ the wall normal |
| $L, t$ | wall length and thickness |
| $\mathbf A, \mathbf B, \mathbf D$ | laminate extensional, coupling and bending stiffness (classical lamination theory) |
| $S_{eq}$ | equivalent transverse shear stiffness of a wall (N/mm) |
| $D_b$ | bending stiffness used by the wall beam (depends on the bending option) |
| $\beta$ | wall shear parameter $12 D_b /(L^2 S_{eq})$ |
| $\psi$ | bend–twist ratio $D_{16}/\sqrt{D_{11} D_{66}}$ of a wall laminate |
| $\lambda$ | twist length $\sqrt{D_{66}/A_{44}}\,/L$ |
| $\eta$ | twist release; $\Pi = \eta\psi^2$ the bend–twist index |
| $\chi_B$ | membrane–bending measure $B_{11}^2/(A_{11}D_{11})$ |
| $\mathbf C^H, \mathbf S$ | homogenised 6×6 stiffness and compliance |

**Voigt order.** Strains and stresses use the engineering order
$\boldsymbol\varepsilon = [\varepsilon_{11}, \varepsilon_{22}, \varepsilon_{33}, \gamma_{23}, \gamma_{13}, \gamma_{12}]^T$,
$\boldsymbol\sigma^H = \mathbf C^H \boldsymbol\varepsilon$.

**Shell stiffness of the source paper (Eq. 5).** The sandwich-core shell matrix uses
$Q_{11}, Q_{22}, Q_{12}, Q_{66}$ (plane stress, $\sigma_{33}=0$) and the transverse shear terms
$Q_{44} = C_{1313} = G_{13}$ and $Q_{55} = C_{2323} = G_{23}$. RVE Studio reports both the full 3D $\mathbf C^H$ and this shell
form; never compare a shell term (e.g. $Q_{11}$, $\sigma_{33}=0$) with a 3D term ($C_{1111}$, $\varepsilon_{33}=0$).

**Units.** Any consistent set. The supplied libraries use MPa and mm (forces in N, stiffness per length in N/mm);
density is kg/m³ and the mass calculator converts the model length unit (mm, cm, m).

## 2. Ply and laminate

For an orthotropic ply in plane stress, with $\nu_{21} = \nu_{12}E_2/E_1$,

$$
Q_{11} = \frac{E_1}{1-\nu_{12}\nu_{21}},\quad Q_{22} = \frac{E_2}{1-\nu_{12}\nu_{21}},\quad
Q_{12} = \frac{\nu_{12}E_2}{1-\nu_{12}\nu_{21}},\quad Q_{66} = G_{12},
$$

rotated to the ply angle $\theta$ (measured from the wall axis $\bar 1$) to give $\bar{\mathbf Q}(\theta)$ in the
engineering-shear convention. For ply interfaces $z_0 < z_1 < \dots < z_n$ through the wall thickness,

$$
\mathbf A = \sum_k \bar{\mathbf Q}_k (z_k - z_{k-1}),\quad
\mathbf B = \tfrac12 \sum_k \bar{\mathbf Q}_k (z_k^2 - z_{k-1}^2),\quad
\mathbf D = \tfrac13 \sum_k \bar{\mathbf Q}_k (z_k^3 - z_{k-1}^3),
$$

and the transverse shear stiffness $A_{44} = \sum_k G_{th,k}(z_k - z_{k-1})$ with $G_{th}$ the through-thickness
shear modulus of each ply (default $G_{12}$).

**Neutral plane.** For an unsymmetric laminate the wall bends about its neutral plane. With
$\mathbf C_{ABD} = \begin{bmatrix}\mathbf A & \mathbf B\\ \mathbf B & \mathbf D\end{bmatrix}^{-1}$ the offset is
$d = -C_{ABD,14}/C_{ABD,44}$ and

$$
\mathbf B_n = \mathbf B - d\,\mathbf A, \qquad \mathbf D_n = \mathbf D - 2d\,\mathbf B + d^2 \mathbf A .
$$

**Lamination parameters and invariants.** The Laminate Tools page reports the five material invariants of each
constituent and the A- and D-lamination parameters $\xi^{A,D}_{1\ldots4}$ on the basis
$[\cos 2\theta, \cos 4\theta, \sin 2\theta, \sin 4\theta]$; hybrid stacks are reported per material.

## 3. Equivalent transverse shear stiffness of a wall

The wall's transverse shear flexibility follows from the shear-stress distribution of a laminated beam. With the
first moment of the axial stiffness about the neutral axis

$$
g(z) = \int_{z_0}^{z} \bar Q_{11}(s)\,(s-d)\,\mathrm ds ,
$$

the equivalent shear stiffness and shear correction factor are

$$
S_{eq} = \frac{D_{n,11}^2}{\displaystyle\int_{-t/2}^{t/2} \frac{g(z)^2}{G_{th}(z)}\,\mathrm dz},
\qquad \kappa = \frac{S_{eq}}{A_{44}} .
$$

$g$ is quadratic in $z$ inside each ply, so the integrand is a quartic; RVE Studio integrates it exactly with
three-point Gauss–Legendre quadrature per ply ($\kappa = 5/6$ for a homogeneous wall).

## 4. The strain-energy wall element

Each wall is a beam of length $L$ and height $h$ with eight degrees of freedom in its local frame,

$$
\mathbf d_e = [\bar u_{1i}, \bar u_{2i}, \bar u_{3i}, \varphi_i, \bar u_{1j}, \bar u_{2j}, \bar u_{3j}, \varphi_j]^T ,
$$

where $\bar u_1$ is along the wall, $\bar u_2$ normal to it in the cell plane, $\bar u_3 = u_3$ along the core
height and $\varphi$ the rotation about $x_3$. Its strain energy, including the homogeneous through-height strain
$\varepsilon_{33}$, is

$$
U_e = \tfrac12\,\mathbf d_e^T \mathbf K_e \mathbf d_e + \mathbf d_e^T \mathbf g_e\,\varepsilon_{33} + \tfrac12 r_e\,\varepsilon_{33}^2 .
$$

* **Membrane part** (wall stretching $\bar\varepsilon_{11}$ and in-plane shear $\bar\gamma_{13}$ of the wall):
  $\dfrac{h}{L}\begin{bmatrix}A_{11} & A_{16}\\ A_{16} & A_{66}\end{bmatrix}$ acting on
  $[\bar u_{1j}-\bar u_{1i},\ \bar u_{3j}-\bar u_{3i}]$.
* **$\varepsilon_{33}$ coupling:** $\mathbf g_e = h\,[-A_{12}, 0, -A_{26}, 0, A_{12}, 0, A_{26}, 0]^T$ and
  $r_e = A_{22} L h$.
* **Bending in the cell plane** (deflection $\bar u_2$, rotation $\varphi$): a two-node beam with bending
  stiffness $D_b h$ and shear stiffness $S_{eq} h$. With $\beta = 12 D_b /(L^2 S_{eq})$ (Timoshenko) or
  $\beta = 0$ (Euler–Bernoulli),

$$
\mathbf K_b = \frac{D_b h}{L^3 (1+\beta)}
\begin{bmatrix}
12 & 6L & -12 & 6L\\
6L & (4+\beta)L^2 & -6L & (2-\beta)L^2\\
-12 & -6L & 12 & -6L\\
6L & (2-\beta)L^2 & -6L & (4+\beta)L^2
\end{bmatrix}.
$$

The element is rotated to the global frame with the wall angle $\alpha = \operatorname{atan2}(y_J-y_I, x_J-x_I)$ and
multiplied by the wall weight $w_e$ (1 for an interior wall, ½ for a wall shared by two periodic images).

**Wall beam option.** *Euler–Bernoulli* ($\beta = 0$) is the original model. *Timoshenko* adds the wall's
transverse shear flexibility; it lowers mainly $Q_{66}$ ($C_{1212}$, $G_{12}$) and $E_1, E_2$, by roughly
$100\beta\,\%$ of the bending-governed terms.

**Wall direction matters.** Ply angles are measured from $\bar 1$, i.e. from Node I towards Node J. Reversing a
wall mirrors ±θ plies and turns an unsymmetric stack over, so for off-axis or unsymmetric laminates the effective
properties depend on the I → J order (no effect for symmetric 0/90 stacks). The benchmark honeycomb has every
wall pointing into the junction node.

## 5. Wall bending options

The wall beam carries only the bending curvature $\kappa_{11}$ about $x_3$. A plate wall can also curve
transversely ($\kappa_{33}$) and twist ($\kappa_{13}$); the bending option decides what happens to those
curvatures and hence the value of $D_b$. *Cylindrical* uses the neutral-plane bending stiffness $D_{n,11}$
(Section 2). The other options start from the condensed bending matrix
$\mathbf D^{\ast} = \mathbf D - \mathbf B\mathbf A^{-1}\mathbf B$ (equal to $\mathbf D$ for a symmetric wall);
in the table $D_{ij}$ stands for $D^{\ast}_{ij}$.

| Option | $D_b$ | Physical meaning |
|---|---|---|
| **Cylindrical** (original model) | $D_{n,11}$ | transverse curvature and twist suppressed (cylindrical bending) |
| **Anisotropic strip** | $D_{11}$ | cylindrical bending with the membrane–bending coupling condensed; differs from Cylindrical only for unsymmetric walls |
| **Free moment** | $D_{11} - [D_{13}\ D_{16}]\begin{bmatrix}D_{33}&D_{36}\\D_{36}&D_{66}\end{bmatrix}^{-1}[D_{13}\ D_{16}]^T$ | transverse and twisting moments free: a lower bound on $D_b$ |
| **Twist-corrected** (calibrated) | $D_{11}(1-\eta\psi^2)$ | the wall releases a calibrated part of its bend–twist coupling |

### 5.1 Why a twist correction is needed

A wall with bend–twist coupling ($D_{16} \ne 0$, e.g. a [θ/−θ] stack that is not balanced through the
thickness, or any off-axis stack) twists when it bends. In the cellular core the wall is held at both ends by the
junctions and over its height by its neighbours, so it is neither free to twist nor fully prevented from twisting.
Cylindrical bending forbids the twist and over-estimates $D_b$ by the factor $1/(1-\psi^2)$ in the limit of full
release; free-moment bending allows all of it and under-estimates $D_b$. The truth lies between, and how much
twist is released depends on how stiff the wall is in twisting relative to its transverse shear over its length.

### 5.2 Twist length and twist release

The **twist length**

$$
\lambda = \frac{\sqrt{D_{66}/A_{44}}}{L}
$$

is the length over which transverse shear resists twisting, relative to the wall length. Short, thick walls have
large $\lambda$ and release more twist; long, thin walls have small $\lambda$ and behave cylindrically. The twist
release

$$
\boxed{\ \eta = \min\!\left[1,\ (a + b\,\psi^2)\,\lambda^{\,n}\right],\qquad a = 2.0,\ b = 2.15,\ n = 0.73\ }
$$

is the fraction of the full bend–twist reduction $\psi^2$ that the wall releases, so

$$
D_b = D_{11}\,(1 - \eta\,\psi^2) .
$$

The coefficients are absolute (no reference wall or configuration). They were fitted to periodic shell FE of
**320 random walls** — 4 to 8 plies of glass, carbon and boron/epoxy with 12 ply angles, $0.031 \le t/L \le 0.124$,
$0.010 \le \lambda \le 0.100$, $|\psi| \le 0.85$ — each solved with three mesh-converged shell meshes
(2.5, 1.25 and 0.625 mm). Five-fold cross-validation and 48 independent walls confirm the fit. With it every
in-plane constant of every tested wall with $\beta \le 0.10$ is within 2 % of the converged FE, up to
$|\psi| = 0.85$.

### 5.3 Editing the coefficients

The coefficients $(a, b, n)$ can be changed in the Solver tab, on Benchmark Designs and in Generate DB
(`twist_release=(a, b, n)` in the API, `--twist A B N` on the command line). Use this to recalibrate the wall
against your own FE or tests. A design solved with coefficients other than the calibrated ones is rated
*use with care* (Section 9), because the verified limits were established with the calibrated values.

## 6. Periodic RVE and homogenisation

### 6.1 Periodicity

For a periodic node pair (master $m$, slave $s$) separated by $\Delta\mathbf x = n_1\mathbf a_1 + n_2\mathbf a_2$,

$$
\begin{aligned}
u_{1s}-u_{1m} &= \varepsilon_{11}\Delta x_1 + \tfrac12\gamma_{12}\Delta x_2,\\
u_{2s}-u_{2m} &= \varepsilon_{22}\Delta x_2 + \tfrac12\gamma_{12}\Delta x_1,\\
u_{3s}-u_{3m} &= \gamma_{13}\Delta x_1 + \gamma_{23}\Delta x_2,\qquad \varphi_s = \varphi_m .
\end{aligned}
$$

Three translations of one anchor node remove the rigid-body motion. RVE Studio detects periodic pairs from the
fractional lattice coordinates; explicit pairs can be entered for unusual RVEs. All constraints are assembled as
$\mathbf C\mathbf u = \mathbf D\boldsymbol\varepsilon$ with linearly dependent rows removed.

### 6.2 Condensation

After assembly the RVE energy is
$U = \tfrac12\mathbf u^T\mathbf K\mathbf u + \mathbf u^T\mathbf H\boldsymbol\varepsilon + \tfrac12\boldsymbol\varepsilon^T\mathbf R\boldsymbol\varepsilon$.
Stationarity under the periodic constraints gives one KKT system with six right-hand sides (one per unit
macroscopic strain),

$$
\begin{bmatrix}\mathbf K & \mathbf C^T\\ \mathbf C & \mathbf 0\end{bmatrix}
\begin{bmatrix}\mathbf U\\ \boldsymbol\Lambda\end{bmatrix} =
\begin{bmatrix}-\mathbf H\\ \mathbf D\end{bmatrix},
$$

and the homogenised stiffness from the condensed energy,

$$
\mathbf C^H = \frac{1}{V}\left(\mathbf U^T\mathbf K\mathbf U + \mathbf U^T\mathbf H + \mathbf H^T\mathbf U + \mathbf R\right),
\qquad \mathbf S = (\mathbf C^H)^{-1} .
$$

A typical RVE has a few tens of degrees of freedom; RVE Studio solves it with a dense LU (sparse LU above 600 unknowns),
about 1 ms per design.

## 7. Effective properties

From $\mathbf S$: $E_1 = 1/S_{11}$, $E_2 = 1/S_{22}$, $E_3 = 1/S_{33}$, $G_{23} = 1/S_{44}$, $G_{13} = 1/S_{55}$,
$G_{12} = 1/S_{66}$ and $\nu_{ij} = -S_{ji}/S_{ii}$.

The **shell (Q) matrix** is $\mathbf C^H$ condensed to $\sigma_{33} = 0$: the in-plane part from the compliance
sub-matrix of indices (1, 2, 6), and the transverse shear part $Q_{44} = C_{1313}$, $Q_{55} = C_{2323}$ (paper
Eq. 5). Specific stiffness is the property divided by the average density, in GPa/(Mg/m³).

**Which properties depend on which wall terms.**

| Property group | Governed by | Sensitive to bending option / wall beam |
|---|---|---|
| $E_1, E_2, G_{12}, \nu_{12}, \nu_{13}$ ($C_{1111}$, $C_{2222}$, $C_{1212}$, $Q_{11}$, $Q_{22}$, $Q_{66}$ in bending-dominated cells) | wall bending $D_b$ | yes |
| $E_3$, $G_{13}$, $G_{23}$, $\nu_{31}$ | wall membrane $\mathbf A$ | no |
| $C_{1111}$, $C_{2222}$, $C_{1122}$ in stretch-dominated cells (triangular) | wall membrane $A_{11}$ | only through $\chi_B$ |

## 8. Physical properties

For wall $e$ and ply $k$ the solid volume is $V_{ek} = w_e L_e h\,t_{ek}$. The relative density is
$\rho^{\ast} = \sum V_{ek}/V$; with known ply densities the RVE mass is $m = \sum \rho_k V_{ek}$, the average density
$\bar\rho = m/V$, and with costs per kg the material cost per RVE and per unit volume.

## 9. Reliability of the wall model

Every solve is rated **reliable**, **use with care** or **outside the verified range**, separately for the
in-plane and the out-of-plane properties, from three measures evaluated for each wall laminate at its shortest
wall length:

| Measure | Definition | Effect | Reliable | Use with care |
|---|---|---|---|---|
| Bend–twist index (cylindrical Timoshenko wall) | $\Pi = \eta\psi^2$ | in-plane properties over-predicted by ≈ $100\,\Pi$ % | $\Pi \le 0.018$ (≈ 2 %) | $\Pi \le 0.045$ (≈ 5 %) |
| Bend–twist ratio (twist-corrected wall) | $\lvert\psi\rvert$, with $0.010 \le \lambda \le 0.100$ | — | $\lvert\psi\rvert \le 0.85$ | λ outside the calibrated range |
| Wall shear and bend–twist, Euler–Bernoulli wall | $\Pi + \beta$ (cylindrical), $\beta$ (twist-corrected); $\beta = 12D_b/(L^2S_{eq})$ | bending terms ≈ $100(\Pi + \beta)$ % too stiff: the two errors add | $\le 0.02$ | $\le 0.05$ |
| Wall shear, Timoshenko wall | $\beta$ | shear model accuracy | $\beta \le 0.10$ (≈ 2 %) | $\beta \le 0.50$ (≈ 5 %) |
| Membrane–bending coupling | $\chi_B = B_{11}^2/(A_{11}D_{11})$ | $C_{1111}, C_{2222}, C_{1122}$ ($Q_{11}, Q_{12}$) over-predicted by ≈ $100\chi_B$ % | $\chi_B \le 0.02$ | $\chi_B \le 0.05$ |

* The out-of-plane properties ($E_3$, $G_{13}$, $G_{23}$, $\nu_{31}$) depend on none of the bending measures and are
  reliable for every symmetric wall (only $\chi_B$ applies).
* *Anisotropic strip* and *Free moment* are bounds, not verified models; their in-plane properties are always
  rated *outside the verified range*.
* Twist-corrected walls solved with user-defined coefficients are rated at best *use with care*.
* Against the converged shell FE of the 320 calibration walls, every wall rated *reliable* is within 2.0 % in all
  in-plane properties and every wall rated *use with care* within 5.2 %, for each of the four combinations of wall
  beam (Euler–Bernoulli, Timoshenko) and bending option (Cylindrical, Twist-corrected).

The rating, the measures of every wall and plain-language notes appear in the Results tab, in
`summary(result)['reliability']`, in the command-line output and in the `PsiMax`, `BendTwistIndexMax`,
`ChiBMax`, `BetaMax`, `ReliabilityInPlane` and `ReliabilityOutOfPlane` columns of every database row.

## 10. Choosing a wall model — guidance

1. **Reproducing the source paper:** Euler–Bernoulli + Cylindrical (the defaults). All 84 original model
   values are reproduced within 0.3 %.
2. **Design with walls of any stacking sequence:** Timoshenko + Twist-corrected. This is the most accurate
   combination against converged shell FE: within 2 % for $|\psi| \le 0.85$ and $\beta \le 0.10$.
3. **Walls without bend–twist coupling** (cross-ply, balanced-symmetric with small $D_{16}$, $\Pi \le 0.018$):
   Cylindrical and Twist-corrected agree; use Timoshenko whenever $\beta > 0.02$.
4. **Unsymmetric walls:** check $\chi_B$; *Anisotropic strip* removes the coupling from $D_b$ and gives the
   bound with B condensed.
5. **Bounds:** Cylindrical (upper) and Free moment (lower) bracket the wall bending stiffness of any laminate.
6. **When the rating is *use with care* or *outside*:** run *FE Design & Validation* for that design.

## 11. Benchmark designs

The 14 benchmark designs are the regular hexagonal core of the source paper (walls 50 mm, $h = 150$ mm,
five equal plies, E-glass/epoxy and boron/epoxy, 14 stacking sequences), with the six shell stiffnesses
$Q_{11}, Q_{22}, Q_{12}, Q_{66}, Q_{44}, Q_{55}$ each — 84 values. The reference values are those of the
cylindrical Euler–Bernoulli model (paper Tables 3–9). The Benchmark Designs page and `rvestudio benchmark` solve the
designs with any wall beam and bending option:

| Wall model | Within ±1 % of the reference | Largest difference |
|---|---|---|
| Euler–Bernoulli, Cylindrical | 84 / 84 | 0.28 % |
| Timoshenko, Cylindrical | 70 / 84 | $Q_{66}$ 1.2–14.5 % lower |
| Euler–Bernoulli, Twist-corrected | 78 / 84 | $Q_{66}$ up to 26.9 % lower (bend–twist coupled walls) |
| Timoshenko, Twist-corrected | 70 / 84 | $Q_{66}$ 2.1–27.6 % lower |

Differences from the reference with other wall models are physical (wall shear flexibility, twist release), not
errors; their accuracy is judged against shell FE (Sections 12 and 13).

## 12. Shell FE validation

*FE Design & Validation* builds a periodic shell model of the same RVE and compares every term with the
analytical model of the same design. The shell model uses four-node shell elements with the same laminate
(layered section, $G_{13} = G_{23} = G_{th}$, shear correction 5/6) — either the built-in solver (no licence,
seconds) or Abaqus (S4 elements, identical node sets, boundary conditions and energy formulae).

**FE RVE.** By default the smallest rectangular block of 1×1 … 3×3 periodic cells that holds a complete cell
strictly inside, with no wall piece shorter than 5 % of the median wall (for the benchmark honeycomb
173.2 × 150 × 150 mm). It can be chosen as n₁ × n₂ cells or drawn over the replicated design. The PBC results do
not depend on the RVE choice; the DBC results do.

**PBC model (strain energy).** Side faces periodic, top and bottom faces periodic with
$u_3^+ - u_3^- = \varepsilon_{33}h$. Six unit macroscopic strains give the full $\mathbf C$ from the strain energy
(in Abaqus 21 load cases, the six strains and their 15 pairs):

$$
C_{ij} = \frac{W_{ij} - W_i - W_j}{V e^2}\ (i \ne j), \qquad C_{ii} = \frac{2W_i}{V e^2},
$$

and from $\mathbf C$ the Q matrix ($\sigma_{33} = 0$) and the engineering constants. Because the PBC model is
exactly periodic in $x_3$, its results are independent of the RVE height; a thin RVE (e.g. $h = 10$ mm) makes very
fine meshes affordable.

**DBC model.** The side-face conditions of the source study (thesis Table 3.4). With $u_3$ free on top and bottom
it gives the shell terms $Q_{11}\ldots Q_{55}$; with $u_3 = 0$ the 3D terms $C_{1111}, C_{2222}, C_{1122}, C_{1212}$.

**Colour code.** Analytical − FE differences ≤ 1 % green, ≤ 5 % amber, > 5 % red.

## 13. Mesh convergence: is my FE result converged?

An FE stiffness converges towards its exact value as the mesh is refined; terms governed by wall bending
($Q_{66}$, $C_{1212}$, $G_{12}$, $E_1$, $E_2$) converge slowly, and walls with strong bend–twist coupling need finer
meshes. A single mesh cannot tell you its own error. RVE Studio therefore asks for **three meshes with a constant
refinement ratio** $r = h_2/h_1 = h_3/h_2$ (normally 2, e.g. 10, 5, 2.5 mm or 2.5, 1.25, 0.625 mm) and judges every
term with the grid-convergence index of Roache and Celik et al. (2008).

**Observed order and extrapolation.** With values $f_1, f_2, f_3$ on the fine, medium and coarse meshes and
$e_{21} = f_2 - f_1$, $e_{32} = f_3 - f_2$, the observed order $p$ solves

$$
p = \frac{1}{\ln r_{21}}\left|\ln\left|\frac{e_{32}}{e_{21}}\right| + \ln\frac{r_{21}^p - 1}{r_{32}^p - 1}\right| ,
$$

and the Richardson extrapolated value is $f_{ext} = (r^p f_1 - f_2)/(r^p - 1)$.

**Grid-convergence index** (safety factor 1.25):

$$
\mathrm{GCI}_{fine} = \frac{1.25\,|e_{21}/f_1|}{r_{21}^p - 1},\qquad
\mathrm{GCI}_{32} = \frac{1.25\,|e_{32}/f_2|}{r_{32}^p - 1},\qquad
R_{as} = \frac{\mathrm{GCI}_{32}}{r_{21}^p\,\mathrm{GCI}_{fine}} .
$$

$R_{as} = 1$ when all three meshes are in the asymptotic range.

**Verdict for each term** (tolerance 1 %):

| Verdict | Condition | What to do |
|---|---|---|
| **converged** | the meshes differ by less than 0.01 %; or monotone, $0.5 \le p \le 4$, $\lvert R_{as} - 1\rvert \le 0.15$ and $\mathrm{GCI}_{fine} \le 1$ % | quote the extrapolated value (or the finest mesh) |
| **refine (GCI)** | monotone and asymptotic, but $\mathrm{GCI}_{fine} > 1$ % | add a finer mesh (halve the finest size) |
| **refine (not asymptotic)** | implausible $p$ or $\lvert R_{as}-1\rvert > 0.15$ | the coarsest mesh is too coarse: drop it and add a finer one |
| **converged (oscillating)** / **refine (oscillating)** | non-monotone triple; judged by the change between the two finest meshes | refine if the change exceeds 1 % |
| **refine (diverging)** | differences grow with refinement | check the model (element distortion, constraints) |

**Practical rules.**

* As a starting point, put at least 5 elements along the shortest wall on the coarsest mesh, and 10 or more for
  bend–twist coupled walls ($|\psi| > 0.3$); the verdict, not the rule of thumb, decides.
* For PBC validation of the in-plane terms use a thin RVE and meshes such as 2.5/1.25/0.625 mm; out-of-plane
  terms ($E_3$, $G_{13}$, $G_{23}$) converge on any mesh.
* Compare the analytical model with the **extrapolated** FE value of a converged term; an unconverged FE value
  over-states stiffness (shell FE converges from above).
* The report (`FE_STUDY_REPORT.html`) lists $p$, $\mathrm{GCI}_{fine}$, $R_{as}$ and the verdict of every term and
  names the terms that need a finer mesh.

## 14. Wall load recovery

For a macroscopic strain $\boldsymbol\varepsilon$ (or stress, $\boldsymbol\varepsilon = \mathbf S\boldsymbol\sigma$) the nodal
displacements are $\mathbf u = \mathbf U\boldsymbol\varepsilon$. In each wall the membrane strains of the element are
$\bar\varepsilon_{11} = (\bar u_{1j} - \bar u_{1i})/L$, $\bar\varepsilon_{33} = \varepsilon_{33}$,
$\bar\gamma_{13} = (\bar u_{3j} - \bar u_{3i})/L$, and the membrane resultants

$$
[\bar N_{11}, \bar N_{33}, \bar N_{13}]^T = \mathbf A_e\,[\bar\varepsilon_{11}, \bar\varepsilon_{33}, \bar\gamma_{13}]^T .
$$

RVE Studio checks equilibrium, periodic compatibility, macro/micro energy equality, element energy, nodal-force and
$\varepsilon_{33}$-conjugate-force identities for every recovery.

## 15. Local wall buckling

Each wall is screened as a simply supported rectangular laminated plate $a = L$, $b = h$ under its recovered
resultants $(\bar N_{11}, \bar N_{33}, \bar N_{13})$ (compression negative), all scaled by one load factor $\lambda$:

$$
(\mathbf K_b + \lambda\,\mathbf K_g)\,\boldsymbol\phi = \mathbf 0,
\qquad w = \sum_{m}\sum_{n} W_{mn}\sin\frac{m\pi x}{a}\sin\frac{n\pi y}{b} .
$$

$\mathbf K_b$ uses $\mathbf D$ (or $\mathbf D^{\ast} = \mathbf D - \mathbf B\mathbf A^{-1}\mathbf B$ for unsymmetric walls) and keeps
$D_{16}, D_{26}$ and the shear load coupling. The smallest positive $\lambda_{cr}$ is reported with the critical wall,
the dominant mode $(m, n)$ and a truncation check against $(m_{max}-2, n_{max}-2)$. $\lambda_{cr} \ge$ the required
factor passes. This is a screening model: junction rotational restraint, face-sheet restraint, imperfections and
post-buckling are not included.

## 16. Design databases and post-processing

**Generation.** Geometry parameters and target relative density are Latin-hypercube sampled; each wall gets an
independent thickness factor, then all thicknesses are scaled so that
$\rho^{\ast} = \sum w_e L_e t_e / A_{RVE}$ matches the target exactly. Every row carries the wall beam, bending option
and twist coefficients, so a database is solved with the wall model it was generated for.

**Post-processing.** Specific properties, anisotropy ratios, Ashby/scaling plots with power-law fits
$y = C\rho^n$, directional $E(\theta)$ from $E_1, E_2, G_{12}, \nu_{12}$, sensitivity (Pearson, Spearman,
standardised regression, η² for categories), Pareto fronts, weighted and robust ($k\sigma$) multi-objective
screening, manufacturing rules and buckling constraints. These screen the solved database; they are not a
continuous optimiser.

## 17. Limitations

* Linear elasticity, small strains, perfectly bonded walls; junctions are rigid points (no fillet or adhesive).
* The wall is a one-dimensional beam in the cell plane; its accuracy is stated in Section 9. Very thick walls
  ($\beta > 0.5$) and $|\psi| > 0.85$ are outside the verified range.
* The twist release is calibrated for $0.010 \le \lambda \le 0.100$ and symmetric walls; unsymmetric walls are
  rated through $\chi_B$.
* Face sheets are not modelled; the core properties are those of the core alone.
* Buckling is linear screening (Section 15).

## 18. References

1. J. Sriharan, M. Dias, D. Fernando, S. Adhikari (2024). A unified approach for the prediction of the effective
   properties of laminated composite cellular core. *Engineering Structures* 299, 117106.
2. J. Sriharan, M. Dias, S. Adhikari, D. Fernando (2025). A semi-analytical model for predicting the shear buckling of
   laminated composite honeycomb cores in sandwich panels. *Composite Structures* 351, 118629.
3. I.B. Celik, U. Ghia, P.J. Roache, C.J. Freitas, H. Coleman, P.E. Raad (2008). Procedure for estimation and
   reporting of uncertainty due to discretization in CFD applications. *J. Fluids Eng.* 130, 078001.
4. P.J. Roache (1998). *Verification and Validation in Computational Science and Engineering.* Hermosa.
5. R.M. Jones (1999). *Mechanics of Composite Materials*, 2nd ed. Taylor & Francis.
6. L.J. Gibson, M.F. Ashby (1997). *Cellular Solids: Structure and Properties*, 2nd ed. Cambridge University Press.
