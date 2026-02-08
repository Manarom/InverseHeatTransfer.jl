using StaticArrays, Test, BenchmarkTools
include("PolynomialWrappers.jl")
import Polynomials, LegendrePolynomials, SpecialPolynomials
import .PolynomialWrappers as PW

@test Polynomials.Polynomial((1.0,2.0,3.0))(0.5) ≈ PW.StandPolyWrapper((1.0,2.0,3.0))(0.5)
@test  Polynomials.ChebyshevT((1.0,2.0,3.0))(0.5) ≈ PW.ChebPolyWrapper((1.0,2.0,3.0))(0.5)
@test  LegendrePolynomials.Pl(0.5 , 3) ≈ PW.LegPolyWrapper((0.0, 0.0, 0.0, 1.0))(0.5)
b_test = SpecialPolynomials.Bernstein{3}((1.0,2.0,3.0))
b_scaled = PW.ScaledPolynomial(PW.BernsteinSymPolyWrapper((1.0,2.0,3.0)), xmin=0.0, xmax = 1.0)
@test b_test(0.5) ≈ b_scaled(0.5)

#checking derivatives
p_st = Polynomials.Polynomial((1.0, 2.0, 3.0))
p_st_der = Polynomials.derivative(p_st)
b_std = PW.StandPolyWrapper((1.0, 2.0, 3.0))
b_std_der = PW.derivative(b_std)

@benchmark Polynomials.derivative($p_st)
@benchmark PW.derivative($b_std)

b_test = SpecialPolynomials.Bernstein{3}((1.0,2.0,3.0))
derivative(b_test)
PW.derivative(PW.BernsteinPolyWrapper((1.0,2.0,3.0)))