### Optimales Flugzeugrumpfdesign durch Topologieoptimierung (Optimizing aircraft fuselage design using topology optimization)
* Geometry :—
  * Generates a cylindrical mesh with nodes on concentric rings and quadrilateral (quad) elements.
  * Mesh is structured, with nodes defined by radius, angular position, and axial position.
  * Elements are quads connecting nodes in a structured grid, representing fuselage skin panels.
  * Panel area is approximated as total surface area divided by the number of elements.
* Finite Element Method :—
  * Implements a 3D FEM solver for linear static analysis with 3 DOFs per node (x, y, z displacements).
  * Supports mechanical loads (point forces), pressure loads (normal to surface), and thermal loads (uniform thermal expansion).
  * Uses a SIMP (Solid Isotropic Material with Penalization) model with a penalization factor p = 3.
  * Computes a simplified constitutive matrix for plane stress, with adjustments for composites.
  * Assembles a sparse global stiffness matrix using Eigen's triplet-based approach.
  * Applies Dirichlet boundary conditions (fixed ends of the cylinder).
  * Solves the linear system using Eigen's `SimplicialLDLT` solver.
  * Computes outputs like compliance, maximum stress (von Mises based on analytical hoop and longitudinal stresses), buckling factor (simplified panel buckling), fatigue life (S-N curve-based), and maximum deflection.
* Topology Optimization (TO) :—
  * Uses a density-based topology optimization approach (SIMP) with element-wise densities variables.
  * Optimizes for minimum compliance subject to a volume fraction constraint.
  * Applies stress and buckling constraints via penalties in the sensitivity analysis.
  * Updates densities using a bisection method to enforce the volume constraint.
  * Includes a convergence check based on density changes.
  * Adjusts densities aggressively when constraints are violated (e.g., stress or buckling limits).
* Uses Aluminum 2024-T3 as the material.
* Sets up a fuselage with realistic parameters (2 m radius, 10 m length, 4 mm skin thickness, 20 frames, 40 stringers).
* Configures the FEM solver and optimizer with constraints (volume fraction = 0.7, stress ≤ yield strength, buckling factor ≥ 0.6, min thickness = 4 mm).
* Written in C++.
