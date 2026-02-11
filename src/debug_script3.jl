using StaticArrays, Polynomials, Plots
include("PolynomialWrappers.jl")

#V = PolynomialWrappers.VanderMatrix(SVector((1.0,2.3,4.5,6.7,8.9)),PolynomialWrappers.StandPolyWrapper{10,Float64})
a = (1.0,1.0,1.0,2.3,4.5,6.7)
p = PolynomialWrappers.TrigPolyWrapper(a)
p_der = PolynomialWrappers.derivative(p)
x = -1.0:1e-2:1.0
ppp = Polynomials.fit(Polynomial,x,p.(x),20)

plot(x,p.(x))
plot!(x,ppp.(x))

ppp_der = Polynomials.derivative(ppp)
plot(x,ppp_der.(x))
plot!(x,p_der.(x))