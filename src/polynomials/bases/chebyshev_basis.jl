

function eval_scaled_poly(ch::ChebPoly{N,T}, x::S, a, b) where {N,T,S}
        R = promote_type(T, S)
        poly_degree(ch) == -1 && return zero(R)
        poly_degree(ch) == 0 &&  return R(ch.coeffs[1]) 
        
        h = b - a

        ξ = 2 * (x - a) / h - one(T)  # Scale to [-1,1] 

        c0 = R(ch.coeffs[N-1])  # T_{N-2} coeff initially (or 0 if N=2)
        c1 = R(ch.coeffs[N])    # T_{N-1} coeff
    
        @inbounds for k in (N-2):-1:1
            tmp = ch.coeffs[k] - c1
            c1 = c0 + 2ξ * c1
            c0 = tmp
        end
        
    return R(c0 + ξ * c1) 
end
eval_poly(ch::ChebPoly, x) = eval_scaled_poly(ch,x,left_scaler(),right_scaler())
derivative_coefficients(p::ChebPoly)  = derivative_coefficients_scaled(p,left_scaler(),right_scaler())
function derivative_coefficients_scaled(p::ChebPoly{N,T}, a , b) where {N,T}
    s = (b - a)/2.0 # need to scale if the default basis is not -1...1
    b = MVector{N - 1, T}(undef)
    # stensile c'i-1 = c'i+1 + 2(i - 1)*ci i = N,N-1,...,2
    if N ≥ 2
        @inbounds  b[N - 1] = 2.0 * (N - 1) * p.coeffs[N]/s 
    end
    if N ≥ 3
        @inbounds   b[N - 2] = 2.0 * (N - 2) * p.coeffs[N - 1]/s 
    end
    @inbounds   for k in N - 3 : -1 : 1
        b[k] = b[k + 2] +  2.0 * k * p.coeffs[k + 1]/s 
    end
    b[1] /= 2.0
    return b.data
end

(poly::ChebPoly)(x::Number) = eval_poly(poly,x)

function cheb_coefs(vals::AbstractArray{<:Number,N}) where {N}
     # type-I DCT, except for size-1 dimensions where we want identity
    kind = map(n -> n > 1 ? FFTW.REDFT00 : FFTW.DHT, size(vals))
    coefs = FFTW.r2r(vals, kind)

    # renormalize the result to obtain the conventional
    # Chebyshev-polnomial coefficients
    s = size(coefs)
    coefs ./= prod(map(n -> n > 1 ? 2(n-1) : 1, s))
    for dim = 1:N
        if size(coefs, dim) > 1
            coefs[CartesianIndices(ntuple(i -> i == dim ? (2:s[i]-1) : (1:s[i]), Val{N}()))] .*= 2
        end
    end

    return coefs
end

#==#