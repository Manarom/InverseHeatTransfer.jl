



derivative_coefficients(p::StandPoly{N,T}) where {N,T} = ntuple(i -> i * p.coeffs[i + 1], N - 1)
"""
    Function to evaluate polynomials
"""
(p::StandPoly)(x::Number) = eval_poly(p, x)
eval_scaled_poly(p::StandPoly,x,_,_) = evalpoly(x, p.coeffs)
eval_poly(p::StandPoly, x) = evalpoly(x, p.coeffs)