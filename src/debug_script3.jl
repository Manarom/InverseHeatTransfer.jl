using StaticArrays
include("PolynomialWrappers.jl")

V = PolynomialWrappers.VanderMatrix(SVector((1.0,2.3,4.5,6.7,8.9)),PolynomialWrappers.StandPolyWrapper{10,Float64})