
#=function eval_scaled_poly(leg::LegPoly{N,T}, x::S, a, b) where {N,T,S}
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
end=#
function eval_scaled_poly(leg::LegPoly{N,T}, x::S, a, b) where {N,T,S}
    R = promote_type(T, S)
    
    N == 0 && return zero(R)
    N == 1 && return R(leg.coeffs[1])
    
    # Scaling to [-1, 1]
    ξ = R(2 * (x - a) / (b - a) - 1)
    
    b_next2 = zero(R) # b_{k+2}
    b_next1 = zero(R) # b_{k+1}
    
    @inbounds for k in N:-1:2
        n = k - 1 # 
        
        # Clanshaw: b_k = c_k + α_k*ξ*b_{k+1} - β_{k+1}*b_{k+2}
        # Clanshaw-Legendre Polynomials: α_n = (2n+1)/(n+1), β_{n+1} = (n+1)/(n+2)

        α = R(2n + 1) / R(n + 1)
        
        # b_{k+2} -> (n+2)/(n+1) * β_{n+2}
        # Clanshaw:
        # β_val = R(n + 1) / R(n + 2) # β_{k+1}

        b_curr = R(leg.coeffs[k]) + α * ξ * b_next1 - (R(k) / R(k+1)) * b_next2
        
        b_next2 = b_next1
        b_next1 = b_curr
    end
    
    # (P_0 = 1): c_1 + α_0*ξ*b_2 - β_1*b_3
    # n=0: α_0 = 1, β_1 = 1/2
    return R(leg.coeffs[1]) + ξ * b_next1 - R(0.5) * b_next2
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
