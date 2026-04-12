# InverseHeatTransfer

[![Build Status](https://github.com/Manarom/InverseHeatTransfer.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Manarom/InverseHeatTransfer.jl/actions/workflows/CI.yml?query=branch%3Amain)

> [!IMPORTANT]  
> **This package is under active development!** Features and API are subject to change. ❗

**InverseHeatTransfer.jl** is a Julia package designed to solve direct and inverse heat conduction problems. It provides numerical solvers for temperature distribution and robust optimization wrappers for parameter estimation.

## Installation

Since this package and its core dependency are not yet in the General registry, install them directly via URL:

```julia
using Pkg
# Required dependency
Pkg.add(url="https://github.com/Manarom/ScaledPolynomials.jl.git")
# This package
Pkg.add(url="https://github.com/Manarom/InverseHeatTransfer.jl.git")
```

## Mathematical Model

The package solves the 1D heat equation in the following form:

$$
\frac{C}{\lambda} \frac{\partial T}{\partial t} = \frac{\partial^2 T}{\partial x^2} + \frac{\lambda'}{\lambda} \left( \frac{\partial T}{\partial x} \right)^2
$$

With initial conditions:

$$
T(x, 0) = T_i(x)
$$

### Physical Parameters:
- $\lambda$ — Thermal conductivity $[W/(m \cdot K)]$
- $C$ — Volumetric heat capacity $[J/(m^3 \cdot K)]$, where $C = C_p \cdot \rho$
- $C_p$ — Specific heat $[J/(kg \cdot K)]$
- $\rho$ — Density $[kg/m^3]$

## Inverse Problem & Regularization

For ill-posed inverse problems, the package implements advanced estimation and regularization techniques:

*   **Tikhonov Regularization:** Used to ensure solution stability by penalizing high-frequency oscillations in the parameter space.
*   **Weighted Regression:** Supports various weighting strategies for the discrepancy function to handle measurement noise:
    *   **Diagonal Weighting:** Independent weights for each data point.
    *   **Proportional Weighting:** Weights scaled according to the magnitude of the measured values.
    *   **AR(1) Weighting Function:** Accounts for first-order autoregressive correlations in measurement errors.

## Sequential Inverse Solver

The library includes a specialized type for the **sequential solving** of inverse problems. This allows for:
*   Simultaneous estimation of $\lambda$ and $C_p$ across multiple experimental datasets.
*   A **joint discrepancy function** that aggregates residuals from several measurements with different heating regimes.
*   Improved parameter identifiability by leveraging diverse thermal loading scenarios within a single optimization framework.

## Numerical Solvers (Direct Problem)

The direct problem is solved using the following Finite Difference Schemes:
1. **Fully Explicit**
2. **Fully Implicit**
3. **Crank-Nicolson**

## Roadmap (TODO)
- [ ] Integration of **Spectral Methods** via `OrdinaryDiffEq.jl`.
- [ ] Advanced **Sensitivity Analysis** for inverse problems.