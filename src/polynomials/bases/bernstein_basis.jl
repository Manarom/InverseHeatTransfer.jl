"""
    eval_monomial(p::BernsteinSymPoly{D,T}, k::Int, x::T) where {D,T}

Evaluates Bernstein polynomial k'th monomial value for x the index of the monomial
goes from ``0 to D - 1``
"""
function eval_scaled_monomial(p::BernsteinSymPoly{D,S}, k::Int, x::T, a::T, b::T) where {D,S,T}

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
eval_monomial(p::BernsteinPoly{D,S},k::Int,x::T) where {D,T,S}  = eval_scaled_monomial(p, k, x, T(0.0),T(1.0))
"""
    bern_max(::BernsteinPoly{D},k::Int)

Returns Bernstein's monomial (maximum_value, maximum_location) tuple
"""
bern_max(::Type{BernsteinPoly{D,T}},k::Int) where {D,T} = bern_max(BernsteinSymPoly{D,T}, k, T(0.0),T(1.0))
bern_max(::BB, k) where BB <: Union{BernsteinPoly{D,T}, BernsteinSymPoly{D,T}} where {D,T} = bern_max(BB, k)

function bern_max(::Type{BernsteinSymPoly{D,T}}, k::Int, a::T = LEFT_SCALER, b::T = RIGHT_SCALER) where {D,T}
    d = D - 1
    s = b - a
    return (^(k, k) * ^(d, - T(d))*^(d - k, d - k) * binomial(d,k), s * k/d + a)
end

bern_max_locations(p::BernsteinSymPoly{N}) where {N} = [bern_max(p,i)[2] for i in 0 : N-1]
bern_max_values(p::BernsteinSymPoly{N}) where {N} = [bern_max(p,i)[1] for i in 0 : N-1]

function eval_scaled_poly(poly::Union{BernsteinSymPoly{N,S},BernsteinPoly{N,S}},x::T,a::T,b::T) where {N,T,S}
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
# BernsteinPoly   - standard for 
left_scaler(::BernsteinPoly) = 0.0
right_scaler(::BernsteinPoly) = 1.0
(p:: BernsteinSymPoly{N,T})(x::T)  where {N,T} = eval_scaled_poly(p, x, T(left_scaler(p)),T(right_scaler(p)))
(p:: BernsteinPoly{N,T})(x::T) where {N,T} = eval_scaled_poly(p, x, T(left_scaler(p)),T(right_scaler(p)))
derivative_coefficients_scaled(p::BB,a,b) where {BB <: Union{BernsteinSymPoly{N,T},BernsteinPoly{N,T}}} where {N,T} = ntuple(i -> (N - 1)*(p.coeffs[i + 1] - p.coeffs[i])/(b - a), N - 1)

derivative_coefficients(p::BB)  where {BB <: Union{BernsteinSymPoly{N,T},BernsteinPoly{N,T}}} where {N,T} = derivative_coefficients_scaled(p,left_scaler(p), right_scaler(p))