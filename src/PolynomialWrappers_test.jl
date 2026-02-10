using StaticArrays, Test, BenchmarkTools, Plots,LinearAlgebra
include("PolynomialWrappers.jl")
import Polynomials, LegendrePolynomials#, SpecialPolynomials
import .PolynomialWrappers as PW

@test Polynomials.Polynomial((1.0,2.0,3.0))(0.5) ≈ PW.StandPolyWrapper((1.0,2.0,3.0))(0.5)
@test  Polynomials.ChebyshevT((1.0,2.0,3.0))(0.5) ≈ PW.ChebPolyWrapper((1.0,2.0,3.0))(0.5)
@test  LegendrePolynomials.Pl(0.5 , 3) ≈ PW.LegPolyWrapper((0.0, 0.0, 0.0, 1.0))(0.5)
#b_test = SpecialPolynomials.Bernstein{3}((1.0,2.0,3.0))
#b_scaled = PW.ScaledPolynomial(PW.BernsteinSymPolyWrapper((1.0,2.0,3.0)), xmin=0.0, xmax = 1.0)
#@test b_test(0.5) ≈ b_scaled(0.5)

#checking derivatives
p_st = Polynomials.Polynomial((1.0, 2.0, 3.0))
p_st_der = Polynomials.derivative(p_st)
b_std = PW.StandPolyWrapper((1.0, 2.0, 3.0))
b_std_der = PW.derivative(b_std)

@benchmark Polynomials.derivative($p_st)
@benchmark PW.derivative($b_std)

#println("Testing bernstein polynomial derivative (default scaling)")
#b_test = SpecialPolynomials.Bernstein{3}((1.0,2.0,3.0))
#b_der_test = derivative(b_test)
#b_der = PW.derivative(PW.BernsteinSymPolyWrapper((1.0,2.0,3.0)))
# @test norm(PW.coeffs(b_der) .- SpecialPolynomials.coefficients(b_der_tes))≈ 0 atol = 1e-12


println("Testing derivative of default scaled polynomials using Polynomial")

x = collect(range(-1.0,1.0,100))
a = (1.0,22.0,3.0,4.56,3.51)

println("Check unscaled polynomials")
for CheckPolyType in (PW.BernsteinSymPolyWrapper, PW.StandPolyWrapper, PW.ChebPolyWrapper, PW.LegPolyWrapper)
    println("Cheking the derivative for $CheckPolyType")
    pw_poly = CheckPolyType(a)
    pw_poly_vals = pw_poly.(x)
    s_fit = Polynomials.fit(Polynomials.Polynomial, x, pw_poly_vals, PW.poly_degree(pw_poly)) 
    s_der = Polynomials.derivative(s_fit)
    pw_poly_der = PW.derivative(pw_poly)
    @test norm(s_der.(x) .- pw_poly_der.(x)) ≈ 0 atol=1e-12
    @show norm(s_der.(x) .- pw_poly_der.(x))
end


#Scaled polynomial derivatives
x_unscaled = collect(200.0:0.1:2000.0)
a = (1.0,22.0,3.0,4.56,3.51)
(x_scaled,x_min,x_max) = PW.scale_x_to_ξ(x_unscaled)

p_scaled = PW.ScaledPolynomial(PW.ChebPolyWrapper, a,  x_min, x_max)
 y_values = p_scaled.(x_unscaled)
p_scaled_der = PW.derivative(p_scaled)
pp = Polynomials.fit(Polynomials.Polynomial,x_unscaled, y_values, length(a) - 1)
y_der_test = Polynomials.derivative(pp).(x_unscaled)

plot(x_unscaled,y_der_test)
plot!(x_unscaled, p_scaled_der.(x_unscaled))

for CheckPolyType in (PW.BernsteinSymPolyWrapper, PW.StandPolyWrapper, PW.ChebPolyWrapper , PW.LegPolyWrapper)
    println("Cheking the derivative for $CheckPolyType")
    p_scaled = PW.ScaledPolynomial(CheckPolyType, a,  x_min, x_max)
    y_values = p_scaled.(x_unscaled)
    @test norm(y_values .- CheckPolyType(a).(x_scaled)) ≈ 0.0 atol =1e-12
    # 
    pp = Polynomials.fit(Polynomials.Polynomial,x_unscaled, y_values, length(a) - 1)
    y_der_test = Polynomials.derivative(pp).(x_unscaled)

    p_scaled_der = PW.derivative(p_scaled)

    @test norm( p_scaled_der.(x_unscaled) .- y_der_test ) ≈ 0.0 atol = 1e-12

end




#=
plot(x_scaled,y_values)

pp = Polynomials.fit(x_scaled,y_values, 2)

plot!(x_scaled,pp.(x_scaled))

pp_der = Polynomials.derivative(pp)

bern_der_unscaled = PW.derivative(bern_scaled.poly)

dern_der_unscaled = PW.ScaledPolynomial(PW.BernsteinSymPolyWrapper, 2.0*bern_der_unscaled.coeffs./(1800.0),200.0, 2000.0)
plot(x_scaled,pp_der.(x_scaled))
plot!(x_scaled,dern_der_unscaled.(x_scaled))

bern_der_scale = PW.derivative(bern_scaled)
plot!(x_scaled,bern_der_scale.(x_scaled))
=#