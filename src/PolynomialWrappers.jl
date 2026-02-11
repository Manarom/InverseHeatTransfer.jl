module PolynomialWrappers
    using LinearAlgebra,Interpolations,Polynomials,LegendrePolynomials,StaticArrays,RecipesBase
    export BernsteinSymPolyWrapper,StandPolyWrapper, LegPolyWrapper,ChebPolyWrapper, ScaledPolynomial
    const LEFT_SCALER = -1.0
    const RIGHT_SCALER = 1.0
    scalers() = (LEFT_SCALER, RIGHT_SCALER)
    scale_span() = RIGHT_SCALER - LEFT_SCALER
    left_scaler() = LEFT_SCALER
    right_scaler() = RIGHT_SCALER
    abstract type AbstractPolyWrapper{P,V,T} end
    """
    derivative_coefficients(::T) where T <: AbstractPolyWrapper

returns NTuple of polynomial derivatives coefficients 
"""
derivative_coefficients(::T) where T <: AbstractPolyWrapper = throw(error("not implemented on $(T)"))

"""
    derivative_coefficients!(::AbstractVector{T}, ::P) where {T, P <: AbstractPolyWrapper{N,T}} where {N}

Fills in-place the vector polynomial coefficients derivative 

"""
function derivative_coefficients!(a::AbstractVector{T}, poly::P) where {T, P <: AbstractPolyWrapper{N,T}} where {N}
    @assert length(a) == N - 1 "Incorrect vector size"
    copyto!(a,derivative_coefficients(poly))
end

"""
    derivative!(p1::AbstractPolyWrapper,p2::AbstractPolyWrapper)

Fills p1 coefficients from p2 derivative in-place
"""
function derivative!(p1::AbstractPolyWrapper,p2::AbstractPolyWrapper) 
    check_derivative_size_consistency(p1,p2)
    derivative_coefficients!(p1.coeffs, p2)
end
coeffs(p::AbstractPolyWrapper) = getfield(p,:coeffs)

scalers(::AbstractPolyWrapper) =  scalers()
scale_span(::AbstractPolyWrapper) = scale_span() 
left_scaler(::AbstractPolyWrapper) = left_scaler()
right_scaler(::AbstractPolyWrapper) = right_scaler()

    """
    check_derivative_size_consistency(::P, ::Q) where {P <: AbstractPolyWrapper{N,V1,R1}, Q <: AbstractPolyWrapper{M,V2,R2}} where {N, M,V1,V2,R1,R2}

Checks type and size consistency of the first argument to be filled from the second rgument derivative
"""
function check_derivative_size_consistency(::P, ::Q) where {P <: AbstractPolyWrapper{N,V1,R1}, Q <: AbstractPolyWrapper{M,V2,R2}} where {N, M,V1,V2,R1,R2} 
        @assert R1 == R2 "Polynomials must be of the same type"
        @assert N == M - 1 "Inconsistent size of coefficients vector, first argument polynomial should have N - 1 coefficients with respect to the second one"
    end    

const POLY_NAMES_TYPES_DICT = Base.ImmutableDict(
        :trig => :TrigPolyWrapper, 
        :leg => :LegPolyWrapper,
        :stand => :StandPolyWrapper ,
        :chebT => :ChebPolyWrapper,
        :bernstein => :BernsteinPolyWrapper,
        :bernsteinsym => :BernsteinSymPolyWrapper
    )
    for (_poly_name, _PolyType) in  POLY_NAMES_TYPES_DICT
        x = String(_poly_name)
        if _poly_name != :bernsteinsym 
            @eval struct $_PolyType{N,T} <: AbstractPolyWrapper{N,T,Symbol($x)}
                coeffs::MVector{N,T}
                $_PolyType(x::Union{NTuple{N,T},M}) where M <: StaticVector{N,T} where {T,N} = new{N,T}(MVector(x))
            end
        else 
            @eval struct $_PolyType{N,T} <: AbstractPolyWrapper{N,T,Symbol($x)}
                    coeffs::MVector{N,T}
                    binoms::SVector{N,T}
                    function  $_PolyType(coeffs::Union{NTuple{N,T},M}) where M <: StaticVector{N,T} where {N,T} 
                        binoms = SVector{N,T}(binomial(N - 1, i) for i=0 : N - 1)
                        new{N,T}(MVector(coeffs),binoms)
                    end
                end
        end  
        @eval  function $_PolyType(x::AbstractVector{T}) where T 
                    N = length(x)
                    $_PolyType(SVector{N}(x))
                end
        @eval (::Type{$_PolyType{N,T}})(x::StaticVector{N,T})  where {N,T} = $_PolyType(x)
        @eval (::Type{$_PolyType{N,T}})()  where {N,T} = $_PolyType(ntuple(_->zero(T),N))
        @eval function (::Type{$_PolyType{N,T}})(x::AbstractVector{T})  where {N,T}
                @assert N == length(x) "Incorrect vector size"
                @assert isa(T, DataType) && T <: Number "Wrong data type"
                return $_PolyType(x)
            end
        @eval function (::Type{$_PolyType{N}})(x::AbstractVector{T})  where {N,T}
                @assert isa(N,Int) "data parameter {N} must be integer"
                @assert N == length(x) "Incorrect vector size"
                return $_PolyType(x)
        end
        @eval function derivative(p::T) where {T<: $_PolyType{N,Q}} where {N,Q}
                    p  |> derivative_coefficients  |> $_PolyType
              end

    end
    Base.copy(p::AbstractPolyWrapper) = typeof(p)(p.coeffs)
    # derivatives StandardBasis
    derivative_coefficients(p::StandPolyWrapper{N,T}) where {N,T} = ntuple(i -> i * p.coeffs[i + 1], N - 1)

    function derivative_coefficients(p::BernsteinSymPolyWrapper{N,T}) where {N,T}
        ntuple(i -> (N - 1)*(p.coeffs[i + 1] - p.coeffs[i])/scale_span(), N - 1)
    end
    function derivative_coefficients(p::ChebPolyWrapper{N,T}) where {N,T}
        s = scale_span()/2.0 # need to scale if 
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
   function derivative_coefficients(p::LegPolyWrapper{N,T}) where {N,T}
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
    function derivative_coefficients(p::TrigPolyWrapper{N,T}) where {N,T}
            is_ends_with_sin = iseven(N)
            b = is_ends_with_sin ? MVector{N + 1, T}(undef) :  MVector{N, T}(undef)  # if ends with sin, one extra term should be added
           fill!(b,zero(T))
            for k = 2 : N
                degree = k - 1
                n = 1 + floor(degree / 2)
                if isodd(k)  # cos((n-1)πx) → -(n-1)π sin(nπx)  
                    b[k - 1] =  - T(π) * (n-1)/scale_span() * p.coeffs[k] 
                else  # sin(nπx) → +nπ cos(nπx)!            
                    b[k + 1] =   T(π) * n/scale_span() * p.coeffs[k]
                end
            end

            return b.data
        end
    struct ScaledPolynomial{Ptype,T} 
        poly::Ptype
        xmin::T
        xmax::T
        function ScaledPolynomial(::Type{PolyType}, coeffs::Union{NTuple{P,T}, AbstractVector{T}}, xmin::T, xmax::T) where {PolyType <: AbstractPolyWrapper, P,T}
            poly = PolyType(coeffs)
            Ptype = typeof(poly)
            @assert xmin < xmax "xmin must be smaller than xmax"
            new{Ptype,T}(poly,xmin,xmax)
        end
        ScaledPolynomial(p::P;xmin = left_scaler(), xmax = right_scaler()) where P <: AbstractPolyWrapper{N,T} where {N,T} = 
                new{P,T}(p, xmin, xmax)
    end
    function derivative(sp::ScaledPolynomial)
        p_der = derivative(sp.poly)
        p_der.coeffs .*= scale_span_ξ_by_x(sp)
        return ScaledPolynomial(p_der,xmin=left_scaler(sp),xmax = right_scaler(sp))
    end
        """
    (sp::ScaledPolynomial{P,T})(x::T) where {P,T}

When calling ScaledPolynomial on argument, it normalizes its input than calls on normalized 
"""
    (sp::ScaledPolynomial{P,T})(x::T) where {P,T} = scale_x_to_ξ(x, sp) |> sp.poly

    left_scaler(p::ScaledPolynomial) = p.xmin
    right_scaler(p::ScaledPolynomial) = p.xmax
    scalers(p::ScaledPolynomial) = p.xmin,p.xmax
    scale_span(p::ScaledPolynomial) = p.xmax - p.xmin
    coeffs(sp::ScaledPolynomial) = coeffs(sp.poly)
    """
    scale_x_to_ξ(x,x_min,x_max)

Takes x supposing it is scaled from x_min to x_max and scales it 
to default scaling from $(left_scaler()) to $(right_scaler())
"""
scale_x_to_ξ(x, x_min, x_max) = left_scaler() + scale_span() * (x - x_min)/(x_max - x_min) 
scale_ξ_to_x(ξ, x_min, x_max) = x_min +  (x_max - x_min) * (ξ - left_scaler())/scale_span()

scale_span_x_by_ξ(p::ScaledPolynomial) = scale_span(p)/scale_span()
scale_span_ξ_by_x(p::ScaledPolynomial) = scale_span()/scale_span(p)
scale_x_to_ξ(x, p::ScaledPolynomial) = scale_x_to_ξ(x, p.xmin, p.xmax)
scale_ξ_to_x(ξ, p::ScaledPolynomial) = scale_ξ_to_x(ξ, p.xmin, p.xmax)


const AnyPoly  = Union{ScaledPolynomial,AbstractPolyWrapper}
"""
    scale_x_to_ξ(x::AbstractVector)

Makes all elements of vector x to fit in range $(LEFT_SCALER)...$(RIGHT_SCALER)
returns normalized vector , xmin and xmax values
All elements of x must be unique
Makes all elements of vector x to fit in range $(LEFT_SCALER)...$(RIGHT_SCALER)
returns normalized vector , xmin and xmax values
All elements of x must be unique

"""
function scale_x_to_ξ(x::AbstractVector)
    ξ  = copy(x)
    return scale_x_to_ξ!(ξ)
end
function scale_x_to_ξ(x::StaticVector) 
    x_min, x_max = extrema(x)
    s= scale_span()/(x_max - x_min)
    a = left_scaler()
    x_new = @.  s*(x - x_min) + a
    return (x_new , x_min, x_max)   
end
function scale_x_to_ξ!(x)
    x_min, x_max = extrema(x)
    s= scale_span()/(x_max - x_min)
    a = left_scaler()
    @. x  = s*(x - x_min) + a
    return (x , x_min, x_max)    
end
"""
    scale_ξ_to_x(normalized_x::AbstractVector, x_min,x_max)

Creates normal vector from one created with [`scale_x_to_ξ`](@ref)` function 
Assumes `normalized_x` as a vector normalized to $(LEFT_SCALER)...$(RIGHT_SCALER)
preformes the revers operation 
"""
function scale_ξ_to_x(ξ::AbstractVector, x_min, x_max)
    x = copy(ξ)
    return scale_ξ_to_x!(x, x_min, x_max)
end
function scale_ξ_to_x!(ξ, x_min, x_max)
    a = left_scaler()
    s = scale_span()
    @. ξ = x_min + (ξ - a)*(x_max - x_min)/s
    return ξ
end
    """
    refill!(sp::ScaledPolynomial{Ptype,T},new_coeffs::NTuple{N,T}) where {Ptype <: AbstractPolyWrapper{N}} where {N,T}

Fills polynomial coefficients from `new_coeffs` 
"""
function refill!(sp::AnyPoly,new_coeffs::NTuple{N,T})  where {N,T}
        @assert parnumber(sp) == N "Incorrect number of coefficients "
        copyto!(coeffs(sp),new_coeffs)
    end
    """
    refill!(sp::ScaledPolynomial{Ptype,T},new_coeffs::NTuple{N,T}, flag) where {Ptype <: AbstractPolyWrapper{N}} where {N,T}

Fills polynomial coefficients by flag, here flag must be range, bitvector or other types which can be used for 
indexing with the same length as the number of polynomial coefficients 
"""
function refill!(sp::AnyPoly,new_coeffs::NTuple{N,T}, flag)  where {N,T}
        v = @view coeffs(sp)[flag]
        copyto!(v , new_coeffs)
    end    
Base.fill!(p::AnyPoly, v) = fill!(coeffs(p),v)


    const SUPPORTED_POLYNOMIAL_TYPES = Base.ImmutableDict([k=>eval(d) for (k,d) in  POLY_NAMES_TYPES_DICT]...)

    poly_name(::P) where P<: AbstractPolyWrapper{N,T,V} where {N,T,V} = V

    poly_name(::Type{P}) where P<: AbstractPolyWrapper{N,T,V} where {N,T,V} = V

    poly_degree(::AbstractPolyWrapper{N}) where {N} = N - 1

    poly_degree(::Type{P}) where P<:AbstractPolyWrapper{N} where {N} = N - 1

    parnumber(::AbstractPolyWrapper{N,T,V}) where {N,T,V} = N
    parnumber(::ScaledPolynomial{Poly}) where {Poly <: AbstractPolyWrapper{N}} where N = N

"""
    Function to evaluate polynomials
"""
eval_poly(p::StandPolyWrapper,x) = evalpoly(x, p.coeffs)

function eval_poly(ch::ChebPolyWrapper{N,T}, x::S) where {N,T,S}
        R = promote_type(T, S)
        poly_degree(ch) == -1 && return zero(R)
        poly_degree(ch) == 0 &&  return R(ch.coeffs[1]) 
        
        h = scale_span()
        a = left_scaler()

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
function eval_poly(leg::LegPolyWrapper{N,T}, x::S) where {N,T,S}
    R = promote_type(T, S)
    n = N - 1
    n <= 0 && return zero(R)
    n == 0 &&  return R(leg.coeffs[1]) 
    
    h = scale_span() # module default scalers
    (a,) = scalers()
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
function eval_monomial(::TrigPolyWrapper,degree,x) 
     degree != 0 || return 1
     n = 1 + floor(degree/2) 
     return isodd(degree) ? sin(n * pi * x/scale_span()) : cos((n - 1) * pi * x/scale_span())
end
function (poly::Union{StandPolyWrapper,ChebPolyWrapper, LegPolyWrapper})(x::Number)
    return eval_poly(poly,x)
end
"""

"""
function eval_monomial(::BernsteinPolyWrapper{D,T},k::Int,x::Number) where {D,T}  # D - number of polynomial coefficients
    d = D - 1 # polynomial degree
    return binomial(d,k)* ^(one(T) - x, d - k) * x^k
end

"""
    eval_monomial(p::BernsteinSymPolyWrapper{D,T}, k::Int, x::T) where {D,T}

Evaluates Bernstein polynomial k'th monomial value for x the index of the monomial
goes from ``0 to D - 1``
"""
function eval_monomial(p::BernsteinSymPolyWrapper{D,T}, k::Int, x::T) where {D,T}
    k < D && k >= 0 || error("incorrect polynomial degree")
    @inbounds begin
        binom = p.binoms[k + 1] #
        s = scale_span() 
        (a, b) = scalers()
        u = (b - x) / s  # 
        v = (x - a) / s
        
        # Power without ^ : Horner's for log(n) muls
        upow = one(T)
        @simd for _ in 1:(D - 1 - k)
            upow *= u
        end
        vpow = one(T)
        @simd for _ in 1:k
            vpow *= v
        end
        
        return binom * upow * vpow
    end
end

"""
    bern_max(::BernsteinPolyWrapper{D},k::Int)

Returns Bernstein's monomial (maximum_value, maximum_location) tuple
"""
function bern_max(::Type{BernsteinPolyWrapper{D,T}},k::Int) where {D,T}
    d = D - 1
    return (^(k, k) * ^(d, -T(d)) * ^(d - k , d - k) * binomial(d, k), k/d)
end
function bern_max(::Type{BernsteinSymPolyWrapper{D,T}}, k::Int, a::T = LEFT_SCALER, b::T = RIGHT_SCALER) where {D,T}
    d = D - 1
    s = b - a
    return (^(k, k) * ^(d, - T(d))*^(d - k, d - k) * binomial(d,k), s * k/d + a)
end
function bern_max(p::BernsteinSymPolyWrapper{D,T}, k) where {D,T}
    d = D - 1
    s = b - a
    b = p.binoms[k + 1]
    return (^(k, k) * ^(d, - T(d))*^(d - k, d - k) * b , s * k/d + a)

end
bern_max_locations(p::BernsteinSymPolyWrapper{N}) where {N} = [bern_max(p,i)[2] for i in 0 : N-1]
bern_max_values(p::BernsteinSymPolyWrapper{N}) where {N} = [bern_max(p,i)[1] for i in 0 : N-1]
"""
    (poly::Union{TrigPolyWrapper{N,T},BernsteinPolyWrapper{N,T},BernsteinSymPolyWrapper{N,T}})(x::T) where {N,T}

"""
function (poly::Union{TrigPolyWrapper{N,T}, BernsteinPolyWrapper{N,T}, BernsteinSymPolyWrapper{N,T}})(x::T) where {N,T}
    #LegendrePolynomials.Pl(x,l) - computes Legendre polynomial of degree l at point x 
    res = zero(T)
    @inbounds for i ∈ 1 : N
        coeff = poly.coeffs[i]
        res += coeff * eval_monomial(poly, i - 1, x) 
    end
    return res 
    #return sum(ntuple(i -> poly.coeffs[i]*eval_poly(poly,i - 1,x),N))
end

function polyfit(::Type{PV},x::V,y::V, N::Int) where {PV <: AbstractPolyWrapper,V <:AbstractVector{T} } where {T}
    @assert is_in_domain(x) "All values of x must be within $(LEFT_SCALER)...$(RIGHT_SCALER)"
    M = length(x)
    @assert length(y) == M "x and y must be of the same size"
    Vand = Matrix{M,N, T, M*N}(undef)
    p = PV(ntuple(i-> i == 1 ? one(T) : zero(T), N))
    #p = PV([i-> i == 1 ? one(T) : zero(T) for i in 1 : N])
    _fill_vander!(Vand,p,x)
    #=for i = 1 : N
        vi = @view Vand[:,i]
        @. vi = p(x)
    end=#
    return Vand
end
is_in_domain(v) = left_scaler() <= minimum(v) && maximum(v) <= right_scaler()
is_in_domain(p::AnyPoly, v) =  left_scaler(p) <= minimum(v) && maximum(v) <= right_scaler(p)

"""
    VanderMatrix{M <: SMatrix , R <: SMatrix, V <: SVector}
    
This type stores the Vandemonde matrix (fundamental matrix of basis functions),
supports various types of internal polynomials 
Structure VanderMatrix has the following fields:
    v - the matrix itself (each column of this matrix is the value of basis function)
    v_unnorm - version of matrix with unnormalized basis vectors (used for annormalized coefficients of polynomial fitting)
    x_first -  first element of the initial vector 
    x_last  -  the last value of the initial vector
    xi - normalized vector 
    poly_type  - polynomial type name (nothing depends on this name)
"""
struct VanderMatrix{N,CN,T,NxCN,CNxCN,P} #M <: SMatrix , R <: SMatrix, V <: SVector}
    v::SMatrix{N,CN,T,NxCN} # matrix of approximating functions 
    v_unnorm::SMatrix{N,CN,T,NxCN} # unnormalized vandermatrix used to convert fitted parameters if necessarys
    # QR factorization matrices
    Q::SMatrix{N,CN,T,NxCN} 
    R::SMatrix{CN,CN,T,CNxCN}
    x_first::T # first element of the initial array
    x_last::T # normalizing coefficient 
    xi::SVector{N,T} # normalized vector-column 
end# struct spec
"""
    bern_max(::VanderMatrix{N,CN,T,NxCN,CNxCN,P}) where {N,CN,T,NxCN,CNxCN,P<:BernsteinPolyWrapper{CN}}

Returns a vector of maximal values of Bernstein basis polynomial basis for particular VanderMatrix
"""
bern_max_values(::VanderMatrix{N,CN,T,NxCN,CNxCN,P}) where {N,CN,T,NxCN,CNxCN,P<:Union{BernsteinPolyWrapper{CN},BernsteinSymPolyWrapper{CN}}} = [bern_max(P,i)[1] for i in 0:CN-1]
bern_max_locations(::VanderMatrix{N,CN,T,NxCN,CNxCN,P}) where {N,CN,T,NxCN,CNxCN,P<:Union{BernsteinPolyWrapper{CN},BernsteinSymPolyWrapper{CN}}} = [bern_max(P,i)[2] for i in 0:CN-1]

"""
    VanderMatrix(x::StaticArray{Tuple{N},T,1},
                      PolyWrapper::Type{P} # = StandPolyWrapper{CN}
                    ) where {N, T, P <:AbstractPolyWrapper{CN,PN}} where {CN,PN}

Input: 
    x  - vector of independent variables (coordinates)
    Val(CN) - vandermatrix column size (degree of polynomial + 1)
    poly_type - polynomial type name, must be member of SUPPORTED_POLYNOMIAL_TYPES

"""
function VanderMatrix(x::StaticVector{N,T},
                      poly_obj::P # = BernsteinSymPolyWrapper{CN,PN}
                    ) where {N, T, P <: AbstractPolyWrapper{CN,PN}} where {CN,PN}
            # N - number of rows, CN - number of columns
            @assert N >= CN "Degree of polynomial must be less or equal the length og x"
            (_xi,x_first,x_last) = scale_x_to_ξ(x)
            V = Matrix{T}(undef,N,CN) 
            Vunnorm = Matrix{T}(undef,N,CN)  # T{length(x)}(MVector{length(x)}(x))
            fill!(poly_obj, zero(PN))
            _fill_vander!(V, poly_obj,_xi)
            _fill_vander!(Vunnorm, ScaledPolynomial(poly_obj, xmin = x_first, xmax = x_last),x)
            NxCN =  N * CN     
            MatrixType  = SMatrix{N, CN, T,NxCN}
            CNxCN =  CN * CN
            RMatrixType = SMatrix{CN, CN, T,CNxCN}
            VectorType = SVector{N,T}
            _V = MatrixType(V)
            (Q,R) = qr(_V)
            VanderMatrix{N,CN,T,NxCN,CNxCN,P}(_V,# Vandermonde matrix
                MatrixType(Vunnorm), #unnormalized vandermatrix
                MatrixType(Q),
                RMatrixType(R),
                x_first, # first element of the initial array
                x_last, # normalizing coefficient 
                VectorType(_xi),  
            )
end
poly_name(::VanderMatrix{N,CN,T,NxCN,CNxCN,P}) where {N,CN,T,NxCN,CNxCN,P} = poly_name(P)
"""
    _fill_vander!(V, poly_obj::AbstractPolyWrapper,xi)

Function to fill the matrix V columns from polynomial basis functions constructor
with argument vector xi 
"""
function _fill_vander!(V, poly_obj::Union{AbstractPolyWrapper{N,T}, ScaledPolynomial{Ptype}},xi::AbstractVector{D}) where {N, T, D, Ptype <: AbstractPolyWrapper{N} }
    @assert size(V,2) == N "wrong size"
    @assert size(V,1) == length(xi) "wrong size"
    VW = @views eachcol(V)
    fill!(poly_obj,zero(D))
    @inbounds for (i,col) ∈ enumerate(VW)
        coeffs(poly_obj)[i] = one(D)                             
        @. col = poly_obj(xi)
        coeffs(poly_obj)[i] = zero(D)
    end  
    return V
end

"""
    is_the_same_x(v::VanderMatrix,x::AbstractVector)

Checks if input `x` is the same as the one used for VanderMatrix creation
"""
function is_the_same_x(vander::VanderMatrix{N,CN,T},x::AbstractVector{T}) where {N,CN,T}
    (length(x) == N && issorted(x) ) || return false
    x_f = first(x)
    vander.x_first == x_f || return false
    x_l = last(x)
    vander.x_last == x_l || return false
    for i in 1:N 
        x[i] == scale_ξ_to_x(vander.xi[i], x_f, x_l)  || return false
    end
    return true
end
"""
    *(V::VanderMatrix,a::AbstractVector)

VanderMatrix object can be directly multiplyed by a vector
"""
Base.:*(V::VanderMatrix,a::AbstractVector) =V.v*a 
"""
    polyfit(V::VanderMatrix{N,CN,T},x::VT,y::VT) where {N,CN,T<:Number,VT<:Vector{T}}

Fits data x - coordinates, y - values using the VanderMatrix
basis function (this coefficients for normalized x-vector)
```julia
    (a,y_fitted, gf) = polyfit(V,x,y)
    y_fitted =  V.v*a # V.v - is the vandermonde matrix
    # for normalized x, which has all values within [-1,1] range
    # if a = [a₁ , ..., aₙ], 
    # e.g if V is for standard basis:
    (xnorm,) = scale_x_to_ξ(x) # returns vector 
    # y_fitted = a₁ + a₂xnorm + ... + aₙxnormⁿ⁻¹
```

Input:
    x - coordinates, [Nx0]
    y - values, [Nx0]
returns tuple with vector of polynomial coefficients, values of y_fitted at x points
and the norm of goodness of fit     
"""
function polyfit(V::VanderMatrix{N,CN,T},x::VT,y::VT) where {N,CN,T<:Number,VT<:Vector{T}}
    yi =  !is_the_same_x(V,x) ? linear_interpolation(x,y)(scale_ξ_to_x(V)) : y
    a =SVector{CN,T}(V.R\(transpose(V.Q)*yi)) # calculating pseudo-inverse
    y_fit = V*a
    goodness_fit = norm(yi .- y_fit)
    return  (a, y_fit, goodness_fit) 
end
"""
    polyfitn(V::VanderMatrix{N,CN,T},x::VT,y::VT) where {N,CN,T<:Number,VT<:Vector{T}}

Fits data x - coordinates, y - values using the VanderMatrix
basis function (coefficients for unnormalized x-vector)

Input:
    x - coordinates, [Nx0]
    y - values, [Nx0]
returns tuple with vector of polynomial coefficients, values of y_fitted at x points
and the norm of goodness of fit  
"""
function polyfitn(V::VanderMatrix{N,CN,T},x::VT,y::VT) where {N,CN,T<:Number,VT<:Vector{T}}
    yi =  !is_the_same_x(V,x) ? linear_interpolation(x,y)(scale_ξ_to_x(V)) : y
    (Q,R) = qr(V.v_unnorm)
    a =SVector{CN,T}(R\transpose(Q)*yi)# calculating pseudo-inverse
    y_fit = V.v_unnorm * a
    goodness_fit = norm(yi .- y_fit)
    return  (a, y_fit, goodness_fit) 
end




scale_ξ_to_x(V::VanderMatrix) = Vector(scale_ξ_to_x.(V.xi, V.x_first,V.x_last))
function triplicate_columns(a::AbstractVector,T)
    return T(repeat(a,1,3))
end

"""
    fill_box_constraint!(lb,ub,::VanderMatrix{N, CN, T, NxCN, CNxCN, P},
                val_bounds::NTuple{2,T}) where {N, CN, T, NxCN, CNxCN, P<:BernsteinSymPolyWrapper}

Evaluates box-boundaries for polynomial coefficients for `BernsteinSymPolyWrapper` 
polynomial basis
"""
function fill_box_constraint!(lb,ub,::VanderMatrix{N, CN, T, NxCN, CNxCN, P},
                val_bounds::NTuple{2,T}) where {N, CN, T, NxCN, CNxCN, P<:BernsteinSymPolyWrapper}
    fill!(lb,first(val_bounds))
    fill!(ub,last(val_bounds))
end

@recipe function f(m::Union{AbstractPolyWrapper,ScaledPolynomial})
    minorgrid--> true
    gridlinewidth-->2
    dpi-->600
    return t->m.(t)
end


#p = plot(title = "Monomials",legend=:top,legend_columns=2, background_color_legend=RGBA(1, 1, 1, 0.0),foreground_color_legend=nothing)
@recipe function f(V::VanderMatrix{N,CN,T,NxCN,CNxCN,P}; infill = true) where {N,CN,T,NxCN,CNxCN,P}
    for (i,c) in enumerate(eachcol(V.v))
        @series begin 
            label:="$(i)"    
            linewidth:=2
            legend := :top
            legend_columns :=2

            foreground_color_legend :=nothing
            if infill
                fillrange:=0
                fillalpha:=0.3
            end
            markershape:=:none
            (V.xi, c)
        end
    end
end
end