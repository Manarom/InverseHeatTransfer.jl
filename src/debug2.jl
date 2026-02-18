include("InverseHeatTransfer.jl")
using Plots, BenchmarkTools
#using .InverseHeatTransfer
import .InverseHeatTransfer as IT
p = IT.ScaledPolynomial(IT.BernsteinSymPoly((1.0,2.3,5.6)), xmin = 200.0, xmax = 1200.0)
lb = (2.0,5.0,45.0)
flag = (true, true, false)
o = IT.OptimizableVariable(p, lb = lb, flag = flag)

x_new = (23,4.5)

IT.refresh!(o, x_new)
x = collect(200.0:1.0:1200.0)
plot(x,o.(x))

IT.count_lower_bound_violations(o)
@benchmark IT.count_lower_bound_violations($o)

@code_warntype IT.count_lower_bound_violations(o)



