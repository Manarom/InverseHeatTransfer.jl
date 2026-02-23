include("InverseHeatTransfer.jl")
using Plots, BenchmarkTools, Test
#using .InverseHeatTransfer
p = InverseHeatTransfer.ScaledPolynomial(InverseHeatTransfer.BernsteinSymPoly((1.0,6.3,5.6)), xmin = 200.0, xmax = 1200.0)
lb = (2.0,5.0,4.0)
ub = (3.0,7.0,5.0)
flag = (true, true, false)
o = InverseHeatTransfer.OptimizableVariable(p, lb = lb, ub = ub)
I = InverseHeatTransfer
@test I.parnumber(o) == 3
@test I.isoptimizable(o)
I.change_flag(o, new_flag = [false, true, false])
@test all(o.flag .|  .![false, true, false])

I.refill!(o, (1.0,6.3,5.6))
@test  all(I.coeffs(o) .== (1.0,6.3,5.6))
I.modify!(o, [1.0,6.3,6.6][o.flag])
@test  I.fview_coeffs(o) == [1.0,6.3,6.6][o.flag]

I.change_flag(o, new_flag = true)
@test I.count_lower_bound_violations(o) == 1
@test I.count_upper_bound_violations(o) == 1
@test I.count_bound_violations(o) == 2

@test all(I.fview_coeffs(o) .== o.p.poly.coeffs)
@test all(I.fview_lb_coeffs(o) .== lb)
@test all(I.fview_ub_coeffs(o) .== ub)
@test all(I.extract_params(o) .== o.p.poly.coeffs)
@test all(I.extract_optimizable_params(o) .== o.p.poly.coeffs)

I.count_lower_bound_violations(o)
@benchmark I.count_lower_bound_violations($o)

@code_warntype I.count_lower_bound_violations(o)



