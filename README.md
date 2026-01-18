# Implementation of 2D Incompressible Navier–Stokes Solver in Python

## Overview

This project implements a **two-dimensional incompressible Navier–Stokes solver** using the **projection (fractional step) method** in Python.  
The solver computes velocity and pressure fields for viscous flow in a domain containing solid obstacles, using finite-difference discretization on a structured Cartesian grid.

The implementation is intended for **educational and academic purposes**, focusing on clarity and numerical correctness rather than computational efficiency.

---

## Governing Equations

The incompressible Navier–Stokes equations solved are:

### Momentum Equations

∂u/∂t + (u · ∇)u = −(1/ρ)∇p + ν∇²u + f

### Continuity Equation

∇ · u = 0

Where:

- **u = (u, v)** is the velocity field  
- **p** is pressure  
- **ρ** is fluid density  
- **ν** is kinematic viscosity  
- **f** is an external body force  

---

## Assumptions

- Incompressible flow  
- Isothermal conditions  
- No buoyancy effects  
- Constant fluid properties  
- Absolute pressure is irrelevant (only pressure gradients are considered)

---

## Numerical Method

### Spatial Discretization

- Uniform Cartesian grid  
- Second-order central difference scheme for:
  - Pressure gradients
  - Divergence
  - Laplacian (viscous diffusion)

### Time Integration

- Explicit forward Euler method

### Pressure–Velocity Coupling

The **projection method** is used:

1. Predict intermediate velocities without pressure
2. Solve the pressure Poisson equation
3. Correct velocities to enforce incompressibility

### Pressure Solver

- Iterative **Gauss–Seidel method** for the Poisson equation

---

## Domain and Geometry Handling

- The computational domain is defined using a **binary obstacle grid**:
  - `1` → solid wall
  - `0` → fluid region
- A coarse geometry is **upscaled using a Kronecker product** to increase resolution
- **No-slip boundary conditions** are enforced by zeroing velocity inside solid cells

---


