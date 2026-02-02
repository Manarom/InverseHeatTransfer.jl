# InverseHeatTransfer

[![Build Status](https://github.com/Manarom/InverseHeatTransfer.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Manarom/InverseHeatTransfer.jl/actions/workflows/CI.yml?query=branch%3Amain)


Formulates the problem to solve the following equation:
$$
  \frac{C}{λ}Tₜ= Tₓₓ + \frac{λ'}{λ}(Tₓ)² \\
   T(x,0) = Tᵢ(x)
$$
where

$λ$ - thermal conductivity, Kg/m^3 * J/(Kg*K)    

$C$ - thermal capacity,    $C = Cp*ρ$    $Cp$ - specific heat, J/(kg*K), $ρ$ - density, kg/m³ 

$Tₜ = ∂T/∂t$  

$Tₓ = ∂T/∂x$ 

$Tₓₓ = ∂²T/∂x²$ 


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


Finite - difference schemes 

1.) Fully explicit

2.) Fullly implicit

3.) Crank-Nicolson