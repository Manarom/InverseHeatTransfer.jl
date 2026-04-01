# InverseHeatTransfer

[![Build Status](https://github.com/Manarom/InverseHeatTransfer.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Manarom/InverseHeatTransfer.jl/actions/workflows/CI.yml?query=branch%3Amain)


Solves direct problem and wrapps the heat transfer problem for the optimization 

$$
  \frac{C}{\lambda}Tₜ= Tₓₓ + \frac{\lambda'}{\lambda}(T_x)² \\
   T(x,0) = T_i(x)
$$

where

$$ 
    \lambda -\ thermal\ conductivity, Kg/m^3 * J/(Kg*K)    \\
    C - thermal \ capacity,    \\
    C = C_p \cdot \rho     \\
    C_p -\ specific \ heat, J/(kg*K), \\
    \rho -\ density,\ kg/m^2 \\
    T_t = \frac{\partial T} {\partial t}  \\
    T_x = \frac{\partial T} {\partial x} \\
    T_{xx} = \frac{\partial ^2 T} {\partial x ^2}\\
 $$ 


Dirichlet conditions:

$$
    T(0,t) = f(t)
    T(H,t) = g(t)
$$

Dirichlet conditions:

$$
    Tₓ(0,t) = f(t)
    Tₓ(H,t) = g(t)
$$

Robin conditions:

$$
    Tₓ(0,t) = f(T)
    Tₓ(H,t) = g(T)
$$


Direct problem Finite - difference schemes 

1.) Fully explicit

2.) Fullly implicit

3.) Crank-Nicolson

TODO:


 - SpectralMethods + OrdinaryDiffEq
 - Sensitivity analysis