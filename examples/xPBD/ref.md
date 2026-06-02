# High-Order xPBD Solver Implementation Guide

**Target Architecture:** Performance-portable Exascale C++ (MATAR / Kokkos)
**Application:** Dynamic relaxation preconditioner for quasistatic Newton-CG solvers (e.g., within the Fierro ecosystem).

## 1. Mathematical Foundations & Requirements

Extending Extended Position-Based Dynamics (xPBD) to arbitrary-order hexahedral elements requires moving away from the spatially constant deformation gradients of linear tetrahedra and adopting an iso-parametric, quadrature-based approach.

- **Iso-parametric Constraints:** Constraints are derived directly from the continuous elastic energy potential (e.g., Neo-Hookean) evaluated via numerical quadrature, separating the geometric discretization from the material's constitutive law.
- **The Mass Lumping Requirement:** xPBD fundamentally requires a diagonal, lumped mass matrix to avoid global linear solves during constraint projection.
- **The High-Order Mass Trap:** Standard exact integration for mass lumping produces zero or negative masses for polynomial degrees $p \ge 2$. You **must** align your nodal basis with **Gauss-Lobatto-Legendre (GLL)** quadrature points to guarantee strictly positive inverse masses ($w_i > 0$).
- **Under-Integration Advantage:** The exact geometric order of the Hessian integrand is $2p$. A GLL rule with $p+1$ points exactly integrates polynomials up to degree $2p-1$. This slight under-integration acts as an implicit, physically consistent form of anti-locking for nearly incompressible materials.

## 2. Boundary Conditions

### Dirichlet (Displacement) Boundaries

Enforced implicitly inside the xPBD loop via **Inverse Mass Masking**.

- During initialization, set the inverse mass of fixed nodes to zero: $w_i = 0$.
- For prescribed moving boundaries, manually update the coordinate $x_i(t)$ prior to the xPBD loop. The constraint projection will pull interior free nodes to comply, treating the $w_i = 0$ nodes as infinitely heavy objects. This avoids warp-divergent `if` statements in the GPU kernels.

### Neumann (Force) Boundaries

Enforced explicitly *before* the xPBD constraint loop.

- Calculate **consistent nodal forces** by integrating applied surface tractions over the boundary faces using 2D GLL surface quadrature. Do not simply divide total force by node count.
- Apply these forces during the explicit kinematic prediction step:
  $v_i \leftarrow v_i + \Delta t  w_i f_i^{ext}$
  $p_i \leftarrow x_i + \Delta t  v_i$
- *Preconditioning Note:* Ramp boundary forces linearly over several timesteps and apply artificial velocity damping to stabilize the dynamic relaxation.

## 3. The Assembly-Free (Matrix-Free) Paradigm

To achieve Exascale efficiency and avoid memory bandwidth starvation, the solver must operate strictly within thread-local registers.

- **Do not allocate or store** the $(p+1)^3 \times 3$ local gradient vectors $\nabla C_j$.
- **Do not allocate or store** the $9 \times 9$ material stiffness tensors $\mathcal{C}$ at quadrature points.
- **On-the-Fly Contraction:** Evaluate spatial gradients dynamically using tensor-product sum factorization. Calculate the deformation gradient $F$, evaluate the material model, perform the contraction, accumulate the scalar result, and immediately discard the intermediate variables.

## 4. MATAR Algorithmic Architecture

xPBD is a Gauss-Seidel solver. To execute safely across highly concurrent threads without atomic race conditions on shared nodes, elements must be grouped using **Graph Coloring** (the dual graph). 

To prevent catastrophic register spilling when processing high-order elements, the xPBD projection loop is divided into a **Split-Dispatch** architecture. 

### Data Structures

```cpp
CArrayDual<double> x;       // Current positions (num_nodes, 3)
CArrayDual<double> p;       // Predicted unconstrained positions (num_nodes, 3)
CArrayDual<double> v;       // Velocities (num_nodes, 3)
CArrayDevice<double> w;     // Inverse lumped masses (num_nodes)
CArrayDevice<double> f_ext; // Consistent nodal forces (num_nodes, 3)

CArrayDevice<double> lambda; // Element multipliers (num_elements)

// Graph Coloring: Row = color, Length = elements in that color
RaggedCArrayDevice<int> element_colors;
```

