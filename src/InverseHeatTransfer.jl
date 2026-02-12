module InverseHeatTransfer
using LinearAlgebra,Reexport
export HeatTransferProblem
# Write your package code here.
include(joinpath(".","solvers", "OneDHeatTransfer.jl"))
@reexport using .OneDHeatTransfer
include(joinpath(".","polynomials", "PolynomialWrappers.jl"))
@reexport using .PolynomialWrappers
end
