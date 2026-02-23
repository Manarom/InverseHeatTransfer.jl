
function eval_scaled_poly(leg::LegPoly{N,T}, x::S, a, b) where {N,T,S}
    R = promote_type(T, S)
    n = N - 1
    (n < 0 || b <= a) && return zero(R)
    n == 0 &&  return R(leg.coeffs[1]) 
    
    h = b - a # module default scalers
    x_norm = 2 * (x - a) / h - one(T)  # Scale to [-1,1] 
    itr = LegendrePolynomials.LegendrePolynomialIterator(x_norm) 
    # nice feature of LegendrePolynomials  - iterator over monomials
    (s,state) = iterate(itr)
    s *= leg.coeffs[1]
    @inbounds for k in 1 : n
        (v,state) = iterate(itr,state)
        s += R(leg.coeffs[k + 1]) * v
    end
    return s
end
eval_poly(p::LegPoly, x) = eval_scaled_poly(p, x, left_scaler(),right_scaler())

function derivative_coefficients(p::LegPoly{N,T}) where {N,T}
    s =  scale_span()/2.0
    b = zeros(MVector{N - 1, T})
    if N ≥ 2
        @inbounds  b[N - 1] =  (2N - 3) * p.coeffs[N] / s # 0.5 * (N * (N - 1)) *
    end
    if N ≥ 3
        @inbounds  b[N - 2] =  (2N - 5) * p.coeffs[N-1] / s #  0.5 * ((N - 1) * (N - 2)) 
    end
    @inbounds  for k = N - 3 : -1 : 1
        b[k] =  (2k - 1) * (b[k + 2] / (2k + 3) + p.coeffs[k + 1]/s)
    end

    return b.data
end
(poly::LegPoly{N,T})(x) where {N,T} = eval_poly(poly,x)
