module PolynomialWrappers
    using LinearAlgebra,Interpolations,Polynomials,LegendrePolynomials,StaticArrays,RecipesBase, FFTW
    import Base.Broadcast: broadcastable
    export BernsteinSymPoly,StandPoly, LegPoly,ChebPoly, ScaledPolynomial
    
    
    const LEFT_SCALER = -1.0
    const RIGHT_SCALER = 1.0
    scalers() = (LEFT_SCALER, RIGHT_SCALER)
    scale_span() = RIGHT_SCALER - LEFT_SCALER
    left_scaler() = LEFT_SCALER
    right_scaler() = RIGHT_SCALER
    abstract type AbstractPoly{P,V,T} end
    """
    derivative_coefficients(::T) where T <: AbstractPoly

returns NTuple of polynomial derivatives coefficients 
"""
derivative_coefficients(::T) where T <: AbstractPoly = throw(error("not implemented on $(T)"))

"""
    derivative_coefficients!(::AbstractVector{T}, ::P) where {T, P <: AbstractPoly{N,T}} where {N}

Fills in-place the vector polynomial coefficients derivative 

"""
function derivative_coefficients!(a::AbstractVector{T}, poly::P) where {T, P <: AbstractPoly{N,T}} where {N}
    copyto!(a,derivative_coefficients(poly))
end

"""
    derivative!(p1::AbstractPoly,p2::AbstractPoly)

Fills p1 coefficients from p2 derivative in-place
"""
function derivative!(p1::AbstractPoly,p2::AbstractPoly) 
    check_derivative_size_consistency(p1,p2)
    derivative_coefficients!(p1.coeffs, p2)
end
coeffs(p::AbstractPoly) = getfield(p,:coeffs)

scalers(::AbstractPoly) =  scalers()
scale_span(::AbstractPoly) = scale_span() 
left_scaler(::AbstractPoly) = left_scaler()
right_scaler(::AbstractPoly) = right_scaler()

broadcastable(x::AbstractPoly) = Ref(x)
    """
    check_derivative_size_consistency(::P, ::Q) where {P <: AbstractPoly{N,V1,R1}, Q <: AbstractPoly{M,V2,R2}} where {N, M,V1,V2,R1,R2}

Checks type and size consistency of the first argument to be filled from the second rgument derivative
"""
function check_derivative_size_consistency(::P, ::Q) where {P <: AbstractPoly{N,V1,R1}, Q <: AbstractPoly{M,V2,R2}} where {N, M,V1,V2,R1,R2} 
        @assert R1 == R2 "Polynomials must be of the same type"
        @assert N == M - 1 "Inconsistent size of coefficients vector, first argument polynomial should have N - 1 coefficients with respect to the second one"
    end    

const POLY_NAMES_TYPES_DICT = Base.ImmutableDict(
        :trig => :TrigPoly, 
        :leg => :LegPoly,
        :stand => :StandPoly ,
        :chebT => :ChebPoly,
        :bernstein => :BernsteinPoly,
        :bernsteinsym => :BernsteinSymPoly
    )
    for (_poly_name, _PolyType) in  POLY_NAMES_TYPES_DICT
        x = String(_poly_name)
        if _poly_name != :bernsteinsym 
            @eval struct $_PolyType{N,T} <: AbstractPoly{N,T,Symbol($x)}
                coeffs::MVector{N,T}
                $_PolyType(x::Union{NTuple{N,T},M}) where M <: StaticVector{N,T} where {T,N} = new{N,T}(MVector(x))
            end
        else 
            @eval struct $_PolyType{N,T} <: AbstractPoly{N,T,Symbol($x)}
                    coeffs::MVector{N,T}
                    binoms::SVector{N,T}
                    function  $_PolyType(coeffs::Union{NTuple{N,T},M}) where M <: StaticVector{N,T} where {N,T} 
                        binoms = SVector{N,T}(binomial(N - 1, i) for i in 0 : N - 1)
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

    bases_folder = joinpath(".","bases")
    include(joinpath(bases_folder,"bernstein_basis.jl"))
    include(joinpath(bases_folder,"chebyshev_basis.jl"))
    include(joinpath(bases_folder,"trig_basis.jl"))
    include(joinpath(bases_folder,"legendre_basis.jl"))
    include(joinpath(bases_folder,"stand_basis.jl"))
    #include(joinpath(bases_folder,"_basis.jl"))

    Base.copy(p::AbstractPoly) = typeof(p)(p.coeffs)
    # derivatives StandardBasis

    struct ScaledPolynomial{Ptype,T} 
        poly::Ptype
        xmin::T
        xmax::T
        function ScaledPolynomial(::Type{PolyType}, coeffs::Union{NTuple{P,T}, AbstractVector{T}}, xmin::T, xmax::T) where {PolyType <: AbstractPoly, P,T}
            poly = PolyType(coeffs)
            Ptype = typeof(poly)
            @assert xmin < xmax "xmin must be smaller than xmax"
            new{Ptype,T}(poly,xmin,xmax)
        end
        ScaledPolynomial(p::P;xmin = left_scaler(), xmax = right_scaler()) where P <: AbstractPoly{N,T} where {N,T} = 
                new{P,T}(p, xmin, xmax)
    end

            """
        (sp::ScaledPolynomial{P,T})(x::T) where {P,T}

    When calling ScaledPolynomial on argument, it normalizes its input than calls on normalized 
    """
    (sp::ScaledPolynomial)(x) = eval_poly(sp, x)
    eval_poly(p::ScaledPolynomial{P,T}, x::T)  where {P,T} = scale_x_to_ξ(x, sp) |> p.poly
    eval_scaled_poly(p::ScaledPolynomial{P,T}, x::T, _, _) = scale_x_to_ξ(x, sp) |> eval_scaled_poly(p.poly,x,p.xmin,p.xmax)

    left_scaler(p::ScaledPolynomial) = p.xmin
    right_scaler(p::ScaledPolynomial) = p.xmax
    scalers(p::ScaledPolynomial) = p.xmin,p.xmax
    scale_span(p::ScaledPolynomial) = p.xmax - p.xmin
    coeffs(sp::ScaledPolynomial) = coeffs(sp.poly)

    function derivative(sp::ScaledPolynomial)
        p_der = derivative(sp.poly)
        p_der.coeffs .*= scale_span_ξ_by_x(sp)
        return ScaledPolynomial(p_der , xmin = left_scaler(sp) , xmax = right_scaler(sp))
    end
    
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


    const AnyPoly  = Union{ScaledPolynomial,AbstractPoly}
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
        refill!(sp::ScaledPolynomial{Ptype,T},new_coeffs::NTuple{N,T}) where {Ptype <: AbstractPoly{N}} where {N,T}

    Fills polynomial coefficients from `new_coeffs` 
    """
    function refill!(sp::AnyPoly,new_coeffs)
            #@assert parnumber(sp) == N "Incorrect number of coefficients "
            copyto!(coeffs(sp),new_coeffs)
        end
        """
        refill!(sp::ScaledPolynomial{Ptype,T},new_coeffs::NTuple{N,T}, flag) where {Ptype <: AbstractPoly{N}} where {N,T}

    Fills polynomial coefficients by flag, here flag must be range, bitvector or other types which can be used for 
    indexing with the same length as the number of polynomial coefficients 
    """
    function refill!(sp::AnyPoly,new_coeffs, flag)
            v = @view coeffs(sp)[flag]
            copyto!(v , new_coeffs)
        end    
    Base.fill!(p::AnyPoly, v) = fill!(coeffs(p),v)

    const SUPPORTED_POLYNOMIAL_TYPES = Base.ImmutableDict([k=>eval(d) for (k,d) in  POLY_NAMES_TYPES_DICT]...)

    poly_name(::P) where P<: AbstractPoly{N,T,V} where {N,T,V} = V
    poly_name(::Type{P}) where P<: AbstractPoly{N,T,V} where {N,T,V} = V

    poly_degree(::AbstractPoly{N}) where {N} = N - 1
    poly_degree(::Type{P}) where P<:AbstractPoly{N} where {N} = N - 1
    
    parnumber(::AbstractPoly{N,T,V}) where {N,T,V} = N
    parnumber(::ScaledPolynomial{Poly}) where {Poly <: AbstractPoly{N}} where N = N

    function polyfit!(p::Union{AbstractPoly{N,T},ScaledPolynomial{PV}},x::V,y::V) where {PV <:AbstractPoly{N,T}, V <:AbstractVector{D} } where {N,T,D}
        @assert is_in_domain(x) "All values of x must be within range"
        M = length(x)
        @assert length(y) == M "x and y must be of the same size"
        Vand = Matrix{D}(undef, M, N)
        _fill_vander!(Vand,p,x)
        refill!(p,Vand\y)
        return p
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
        bern_max(::VanderMatrix{N,CN,T,NxCN,CNxCN,P}) where {N,CN,T,NxCN,CNxCN,P<:BernsteinPoly{CN}}

    Returns a vector of maximal values of Bernstein basis polynomial basis for particular VanderMatrix
    """
    bern_max_values(::VanderMatrix{N,CN,T,NxCN,CNxCN,P}) where {N,CN,T,NxCN,CNxCN,P<:Union{BernsteinPoly{CN},BernsteinSymPoly{CN}}} = [bern_max(P,i)[1] for i in 0:CN-1]
    bern_max_locations(::VanderMatrix{N,CN,T,NxCN,CNxCN,P}) where {N,CN,T,NxCN,CNxCN,P<:Union{BernsteinPoly{CN},BernsteinSymPoly{CN}}} = [bern_max(P,i)[2] for i in 0:CN-1]

    """
        VanderMatrix(x::StaticArray{Tuple{N},T,1},
                        Poly::Type{P} # = StandPoly{CN}
                        ) where {N, T, P <:AbstractPoly{CN,PN}} where {CN,PN}

    Input: 
        x  - vector of independent variables (coordinates)
        Val(CN) - vandermatrix column size (degree of polynomial + 1)
        poly_type - polynomial type name, must be member of SUPPORTED_POLYNOMIAL_TYPES

    """
    function VanderMatrix(x::StaticVector{N,T},
                        poly_obj::P # = BernsteinSymPoly{CN,PN}
                        ) where {N, T, P <: AbstractPoly{CN,PN}} where {CN,PN}
                # N - number of rows, CN - number of columns
                @assert N >= CN "Degree of polynomial must be less or equal the length og x"
                (_xi,x_first,x_last) = scale_x_to_ξ(x)
                V = Matrix{T}(undef,N,CN) 
                Vunnorm = Matrix{T}(undef,N,CN)  # T{length(x)}(MVector{length(x)}(x))
                fill!(poly_obj, zero(PN))
                _fill_vander!(V, poly_obj,_xi)
                _fill_vander_unscaled!(Vunnorm, poly_obj, x,  x_first, x_last)
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
        _fill_vander!(V, poly_obj::AbstractPoly,xi)

    Function to fill the matrix V columns from polynomial basis functions constructor
    with argument vector xi 
    """
    function _fill_vander!(V, poly_obj::Union{AbstractPoly{N,T}, ScaledPolynomial{Ptype}},xi::AbstractVector{D}) where {N, T, D, Ptype <: AbstractPoly{N} }
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
    function _fill_vander_unscaled!(V , p::AbstractPoly{N,T}, xi::AbstractVector{D} , xmin , xmax) where {N, T, D }
        @assert size(V,2) == N "wrong size"
        @assert size(V,1) == length(xi) "wrong size"
        VW = @views eachcol(V)
        fill!(p, zero(D))
        #f = Base.Fix1(eval_scaled_poly,p)
        @inbounds for (i,col) ∈ enumerate(VW)
            coeffs(p)[i] = one(D)                             
            @. col = eval_scaled_poly(p, xi, xmin, xmax)
            coeffs(p)[i] = zero(D)
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
                    val_bounds::NTuple{2,T}) where {N, CN, T, NxCN, CNxCN, P<:BernsteinSymPoly}

    Evaluates box-boundaries for polynomial coefficients for `BernsteinSymPoly` 
    polynomial basis
    """
    function fill_box_constraint!(lb,ub,::VanderMatrix{N, CN, T, NxCN, CNxCN, P},
                    val_bounds::NTuple{2,T}) where {N, CN, T, NxCN, CNxCN, P<:BernsteinSymPoly}
        fill!(lb,first(val_bounds))
        fill!(ub,last(val_bounds))
    end

    @recipe function f(m::Union{AbstractPoly,ScaledPolynomial})
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