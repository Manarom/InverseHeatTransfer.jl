module InverseHeatTransfer

# Write your package code here.
include(joinpath(".","solvers", "OneDHeatTransfer.jl"))
include(joinpath(".","polynomials", "PolynomialWrappers.jl"))
end
