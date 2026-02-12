function derivative_coefficients(p::BernsteinSymPoly{N,T}) where {N,T}
    ntuple(i -> (N - 1)*(p.coeffs[i + 1] - p.coeffs[i])/scale_span(), N - 1)
end


"""

"""
function eval_monomial(::BernsteinPoly{D,S},k::Int,x::T) where {D,T,S}  # D - number of polynomial coefficients
    R = promote_type(T, S)
    d = D - 1 # polynomial degree
    return binomial(d,k)* ^(one(R) - x, d - k) * x^k
end
"""
    eval_monomial(p::BernsteinSymPoly{D,T}, k::Int, x::T) where {D,T}

Evaluates Bernstein polynomial k'th monomial value for x the index of the monomial
goes from ``0 to D - 1``
"""
function eval_scaled_monomial(p::BernsteinSymPoly{D,S}, k::Int, x::T, a::T, b::T) where {D,S, T}

    @inbounds begin
        R = promote_type(T, S)
        binom = R(p.binoms[k + 1]) #
        s = b - a 
        u = R((b - x) / s)  # 
        v = R((x - a) / s)
       
        upow = one(R)
        @simd for _ in 1:(D - 1 - k)
            upow *= u
        end
        vpow = one(R)
        @simd for _ in 1:k
            vpow *= v
        end
        
        return binom * upow * vpow
    end
end
eval_monomial(p::BernsteinSymPoly{N}, k::Int, x::T) where {N,T} = eval_scaled_monomial(p, k, x, left_scaler(), right_scaler())
"""
    bern_max(::BernsteinPoly{D},k::Int)

Returns Bernstein's monomial (maximum_value, maximum_location) tuple
"""
function bern_max(::Type{BernsteinPoly{D,T}},k::Int) where {D,T}
    d = D - 1
    return (^(k, k) * ^(d, -T(d)) * ^(d - k , d - k) * binomial(d, k), k/d)
end
function bern_max(::Type{BernsteinSymPoly{D,T}}, k::Int, a::T = LEFT_SCALER, b::T = RIGHT_SCALER) where {D,T}
    d = D - 1
    s = b - a
    return (^(k, k) * ^(d, - T(d))*^(d - k, d - k) * binomial(d,k), s * k/d + a)
end
function bern_max(p::BernsteinSymPoly{D,T}, k) where {D,T}
    d = D - 1
    s = b - a
    b = p.binoms[k + 1]
    return (^(k, k) * ^(d, - T(d))*^(d - k, d - k) * b , s * k/d + a)

end
bern_max_locations(p::BernsteinSymPoly{N}) where {N} = [bern_max(p,i)[2] for i in 0 : N-1]

bern_max_values(p::BernsteinSymPoly{N}) where {N} = [bern_max(p,i)[1] for i in 0 : N-1]

eval_poly(p::BernsteinSymPoly{N,S}, x::T) where {N,T,S} = eval_scaled_poly(p,x, left_scaler(),right_scaler())

function eval_scaled_poly(poly::BernsteinSymPoly{N,S},x::T,a::T,b::T) where {N,T,S}
    R = promote_type(S, T)
    t = R( (x - a) / (b - a) )
    beta = MVector{N,R}(poly.coeffs)  # values in this vector are overridden
    @inbounds for j in 1 : (N - 1)
        for k in 1 : (N - j)
            beta[k] = beta[k] * (1 - t) + beta[k + 1] * t
        end
    end
    return beta[1]
end

(poly:: BernsteinSymPoly)(x::Number) = eval_poly(poly,x)