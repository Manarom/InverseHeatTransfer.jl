

function eval_scaled_monomial(::TrigPoly,degree,x,a,b) 
     degree != 0 || return 1
     n = 1 + floor(degree/2) 
     return isodd(degree) ? sin(n * pi * (x - a)/(b - a)) : cos((n - 1) * pi * (x - a)/(b - a))
end

eval_monomial(p::TrigPoly,degree,x) = eval_scaled_monomial(p,degree,x,left_scaler(), right_scaler())
derivative_coefficients(p::TrigPoly) = derivative_coefficients_scaled(p, left_scaler(),right_scaler())
function derivative_coefficients_scaled(p::TrigPoly{N,T},a,b) where {N,T}
    is_ends_with_sin = iseven(N)
    buf = is_ends_with_sin ? MVector{N + 1, T}(undef) :  MVector{N, T}(undef)  # if ends with sin, one extra term should be added
    fill!(buf,zero(T))
    s = b - a
    @inbounds for k = 2 : N
        degree = k - 1
        n = 1 + floor(degree / 2)
        if isodd(k)  # cos((n-1)πx) → -(n-1)π sin(nπx)  
            buf[k - 1] =  - T(π) * (n - 1) * p.coeffs[k] /s
        else  # sin(nπx) → +nπ cos(nπx)!            
            buf[k + 1] =   T(π) * n * p.coeffs[k]/s
        end
    end
    return buf.data
end

"""
    (poly::TrigPoly{N,S})(x::T) where {N,T,S}

"""
function (poly::TrigPoly{N,S})(x::T) where {N,T,S}
    #LegendrePolynomials.Pl(x,l) - computes Legendre polynomial of degree l at point x 
    R = promote_type(S,T)
    res = zero(R)
    @inbounds for i ∈ 1 : N
        res += poly.coeffs[i] * eval_monomial(poly, i - 1, x) 
    end
    return res 
    #return sum(ntuple(i -> poly.coeffs[i]*eval_poly(poly,i - 1,x),N))
end 
function eval_scaled_poly(p::TrigPoly{N,D},x::T, a::T, b::T) where {N,D,T}
    #LegendrePolynomials.Pl(x,l) - computes Legendre polynomial of degree l at point x 
    res = zero(T)
    @inbounds for i ∈ 1 : N
        res += p.coeffs[i] * eval_scaled_monomial(p, i - 1, x, a, b) 
    end
    return res 
end
