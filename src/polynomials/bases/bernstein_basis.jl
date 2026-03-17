"""
    eval_monomial(p::BernsteinSymPoly{D,T}, k::Int, x::T) where {D,T}

Evaluates Bernstein polynomial k'th monomial value for x the index of the monomial
goes from ``0 to D - 1`` where D is the number of coefficients
"""
function eval_scaled_monomial(::BernsteinSymPoly{D,S}, k::Int, x::T, a::T, b::T) where {D,S,T}

    @inbounds begin
        R = promote_type(T, S)
        binom = binomial(D - 1, k) #R(p.binoms[k + 1]) #binomial(N - 1, i)
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

function eval_scaled_poly_old_version(poly::Union{BernsteinSymPoly{N,S},BernsteinPoly{N,S}},x::T,a::T,b::T) where {N,T,S}
    R = promote_type(S, T)
    t = R( (x - a) / (b - a) )
    _one_m_t = one(R) - t
    beta = MVector{N,R}(poly.coeffs)  # values in this vector are overridden
    @inbounds for j in 1 : (N - 1)
        for k in 1 : (N - j)
            beta[k] = beta[k] * _one_m_t  + beta[k + 1] * t
        end
    end
    return beta[1]
end

@inline function eval_scaled_poly(poly::Union{BernsteinSymPoly{N,S},BernsteinPoly{N,S}}, x::T, a::T, b::T) where {N,T,S}
    R = promote_type(S, T)
    t = R((x - a) / (b - a))
    _1_t = one(R) - t
    coeffs = SVector{N, R}(poly.coeffs)
    return _de_casteljau_recursive(coeffs, t, _1_t)
end

# Базовый случай: когда остался один коэффициент — это ответ
@inline _de_casteljau_recursive(coeffs::SVector{1, R}, t, _1_t) where R = coeffs[1]


@inline function _de_casteljau_recursive(coeffs::SVector{K, R}, t, _1_t) where {K, R}
    next_coeffs = SVector{K - 1, R}(ntuple(i -> coeffs[i] * _1_t + coeffs[i+1] * t, Val(K-1)))
    return _de_casteljau_recursive(next_coeffs, t, _1_t)
end

# BernsteinPoly   - standard for 
left_scaler(::BernsteinPoly) = 0.0
right_scaler(::BernsteinPoly) = 1.0
(p:: BernsteinSymPoly{N,T})(x::T)  where {N,T} = eval_scaled_poly(p, x, T(left_scaler(p)),T(right_scaler(p)))
(p:: BernsteinPoly{N,T})(x::T) where {N,T} = eval_scaled_poly(p, x, T(left_scaler(p)),T(right_scaler(p)))
derivative_coefficients_scaled(p::BB,a,b) where {BB <: Union{BernsteinSymPoly{N,T},BernsteinPoly{N,T}}} where {N,T} = ntuple(i -> (N - 1)*(p.coeffs[i + 1] - p.coeffs[i])/(b - a), N - 1)
derivative_coefficients(p::BB)  where {BB <: Union{BernsteinSymPoly{N,T},BernsteinPoly{N,T}}} where {N,T} = derivative_coefficients_scaled(p,left_scaler(p), right_scaler(p))

elevate_degree(p::R) where R <:BernsteinSymPoly{N,T} where {N,T} = p.coeffs |>bernstein_elevate_degree |> BernsteinSymPoly{N + 1, T}  

function bernstein_elevate_degree(b::StaticVector{N,T}) where {N,T}
    b_new = MVector{N + 1,T}(undef)
    n = N - 1
    b_new[1] = b[1]
    for i in 1 : N - 1
        a = i / (n + 1)
        b_new[i + 1] = a *  b[i]  + (1 - a) * b[i + 1] 
    end
    b_new[end] = b[end]
    return b_new
end
function bernstein_elevate_degree!(b, k::Int) 
    for _ in 1 : k
        bernstein_elevate_degree!(b)
    end
    return b 
end 
bernstein_elevate_degree(b::Vector) = copy(b) |> bernstein_elevate_degree!
"""
    bernstein_elevate_degree!(b::Vector)

Elevates the degree of polynomial bernstein bsis by one,
inout argument  - expantion coefficients 
ΣCᵢBᴺᵢ  to ΣC'ᵢBᴺ⁺¹ᵢ 
"""
function bernstein_elevate_degree!(b::AbstractVector)
    N = length(b)
    resize!(b, N + 1)
    n = N - 1
    b[end] = b[end - 1]
    b_i = b[1]
    b_ip1 = b[1]
    @inbounds for i in 1 : N - 1
        a = i / (n + 1)
        b_ip1 = b[i + 1]
        b[i + 1] = a *  b_i  + (1 - a) * b_ip1
        b_i = b_ip1 
    end
    return b
end

#strores Bernstein to standard basis conversion matrix

const BERN_TO_STAND = Dict{Int,Matrix{Float64}}()

function get_bern_to_stand_matrix(N::Int)
    if !haskey(BERN_TO_STAND, N) 
        BERN_TO_STAND[N] = bernstein_to_standard_matrix(N)
    end
    return  BERN_TO_STAND[N]
end

function to_standard_basis(p::BernsteinPoly{N}) where {N}
    return get_bern_to_stand_matrix(N)*p.coeffs
end


"""
    bernstein_matrix(N::Int, T::Type=Float64)

Gives matrix of bernstein to standard basis conversion coefficients
"""
function bernstein_to_standard_matrix(M::Int) 
    mat = zeros(Float64, M, M)
    N = M - 1
    for j in 0 : N
        cnj = binomial(N, j)
        for k in 0 : j
            # T[j,k] = C(N, j) * C(j, k) * (-1)^(j-k)
            val = cnj * binomial(j, k)
            if isodd(j - k)
                mat[j + 1, k + 1] = -val
            else
                mat[j + 1, k + 1] = val
            end
        end
    end
    return mat
end


