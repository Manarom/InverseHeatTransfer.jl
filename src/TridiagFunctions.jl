
"""
    fill_tridiag!(M::Tridiagonal{T,Vector{T}},Fm1::V,F::V,Fp1::V,a0,am1,a,ap1) where V <: AbstractVector{T} where T


Fills three-diagonal matrix from three vectors combination:

upper diagonal = ap1 * Fm1
main diagonal = a0 + a * F
lower diagonal = am1 * Fm1
"""
function fill_tridiag!(M::Tridiagonal{T,Vector{T}},
                                    Fm1::AbstractVector{T} ,F::AbstractVector{T} ,
                                    Fp1::AbstractVector{T} ,a0,am1,a,ap1)  where T

    @. M.dl = am1*Fm1
    @. M.d = a0 + a*F
    @. M.du = ap1*Fp1
    return nothing
end
function fill_tridiag!(M::Tridiagonal{T,Vector{T}},Fm1::AbstractVector{T},F::AbstractVector{T},Fp1::AbstractVector{T},a0,am1,a,ap1,k::Int)  where T
    M.dl[k] = am1*Fm1[k]
    M.d[k] = a0 + a*F[k]
    M.du[k] = ap1*Fp1[k]
    return nothing
end
function fill_tridiag_sym!(M::Tridiagonal{T,Vector{T}},F::AbstractVector{T}, a0,am1,a,ap1, k::Int)  where T
    M.d[k] = a0 + a * F[k]
    M.dl[k] = am1 * F[k]
    M.du[k] = ap1 * F[k]
    return nothing
end

function tridiag_sum(Fm1::AbstractVector{T},F::AbstractVector{T}, Fp1::AbstractVector{T}, a0, am1, a, ap1, k::Int) where T
    return a0 + am1*Fm1[k] +  a*F[k] + ap1*Fp1[k]
end

"""
    tridiag_vect_prod!(b::AbstractVector{T},Fm1::AbstractVector{T},
                                    F::AbstractVector{T}, Fp1::AbstractVector{T}, a0, 
                                    am1, a, ap1,N::Int)  where T

fills the vector b = M*b as M = Tridiagonal(Fm1*am1, a0 + a*F, Fp1*ap1) , N is the length of b

"""
function tridiag_mul!(b::AbstractVector{T},
                                    Fm1::AbstractVector{T},F::AbstractVector{T}, Fp1::AbstractVector{T}, 
                                    a0, am1, a, ap1)  where T
    N = length(F)
    (N == length(b) &&  N - 1 == length(Fp1) && N - 1 == length(Fm1)) || throw(DimensionMismatch("incosistent dimentions"))

    bm1 = b[1]

    b[1] = (a0 + a * F[1]) * b[1] + ap1 * Fp1[1] * b[2]
    @inbounds @fastmath for ii in 2 : N - 1
        b1, b2, b3 = bm1, b[ii], b[ii + 1]
        f1, f2, f3 = Fm1[ii - 1],F[ii],Fp1[ii] 
        b[ii] = am1 * b1 * f1 + b2*(a0 + a*f2) + ap1 *b3 *f3
        bm1 = b2
    end
    b[N] = am1 * bm1 * Fm1[N - 1]  + b[N]*(a0 + a*F[N]) 

    return nothing
end

"""
    column_sym_tridiag_mul!(b::AbstractVector{T}, F::AbstractVector{T}, a0, a, a1)  where T

Fills vector b with the result of its product with square tridiagonal matrix with `column symmetry`

`Column symmetry` means that the matrix has the following structure of `i`'th row: 

``M[i,:] = [..., a1*F[i], a0 + a*F[i], a1*F[i], ...]``

``M[1,:] = [ a0 + a*F[1], a1*F[1], ...]``

``M[end,:] = [..., a1*F[end], a0 + a*F[end]]``

b <= M*b as M = Tridiagonal(F*ap1, a0 + a*F, F*ap1)

!!! modifies b
"""
function column_sym_tridiag_mul!(b::AbstractVector{T}, F::AbstractVector{T}, a0, a, a1)  where T
    N = length(F)
    (N > 0 && N == length(b) ) || throw(DimensionMismatch("incosistent dimentions"))
    bm1 = b[1]
    b[1] = (a0 + a*F[1])*b[1] + a1*F[1]*b[2]

    @inbounds @fastmath  for ii in 2 : N - 1
        b1, b2, b3 = bm1, b[ii], b[ii + 1]
        f = F[ii]
        b[ii] = a1 * (b1 + b3) * f + b2 * ( a0 + a * f )
        bm1 = b2
    end

    b[N] = b[N] * ( a0 + a * F[N] ) + a1 * bm1 * F[N] 

    return nothing
end
"""
    column_sym_tridiag_muladd!(b::AbstractVector{T},c::AbstractVector{T},
                                    F::AbstractVector{T}, a0, a, a1)  where T



Adds to vector b the result of square tridiagonal matrix with `column symmetry` multiplication by vector ```c```

`Column symmetry` means that the matrix has the following structure of `i`'th row: 

``M[i,:] = [..., a1*F[i], a0 + a*F[i], a1*F[i], ...]``

``M[1,:] = [ a0 + a*F[1], a1*F[1], ...]``

``M[end,:] = [..., a1*F[end], a0 + a*F[end]]``

b =b + M*c as M = Tridiagonal(F*ap1, a0 + a*F, F*ap1)

"""
function column_sym_tridiag_muladd!(b::AbstractVector{T},c::AbstractVector{T},
                                    F::AbstractVector{T}, a0, a, a1)  where T
    N = length(F)
    (N > 1 && N == length(c) && N == length(b)) || throw(DimensionMismatch("incosistent dimentions"))

    b[1] += (a0 + a * F[1]) * c[1] + a1 * F[1] * c[2]


    @inbounds @simd for ii in 2 : N - 1   # it is safe to use @simd here since all loops iterations can be run in parallel
        c1, c2, c3 = c[ii - 1], c[ii], c[ii + 1]
        f = F[ii]
        b[ii] += a1*(c1 + c3)*f + c2*(a0 + a*f)
    end

    b[N] += c[N] * ( a0 + a*F[N] ) + a1 * c[N - 1] * F[N] 

    return nothing
end

"""
    tridiag_vect_prod!(b::AbstractVector{T},Fm1::AbstractVector{T},
                                    F::AbstractVector{T}, Fp1::AbstractVector{T}, a0, 
                                    am1, a, ap1)  where T

fills the vector b <= M*c + b as M = Tridiagonal(Fm1*am1, a0 + a*F, Fp1*ap1)
Only `b` is modified
"""
function tridiag_muladd!(b::AbstractVector{T},c::AbstractVector{T}, Fm1::AbstractVector{T},
                                    F::AbstractVector{T}, Fp1::AbstractVector{T}, a0, 
                                    am1, a, ap1)  where T
    N = length(F) 
    (N==length(c) && N == length(b) && N - 1 == length(Fp1) && N - 1 == length(Fm1)) || throw(DimensionMismatch("incosistent dimentions"))


    b[1] += (a0 + a*F[1])*c[1] + ap1 * Fp1[1]*c[2]

    @inbounds @simd for ii in 2 : N - 1
        c1, c2, c3 = c[ii - 1], c[ii], c[ii + 1]
        f1, f2, f3 = Fm1[ii - 1], F[ii], Fp1[ii] 
        b[ii] += am1 * c1 * f1+ c2 * (a0 + a * f2) + ap1 * c3 * f3
    end

    b[N] += am1 * c[N - 1] * Fm1[N - 1]  +  c[N] * (a0 + a * F[N]) 

    return nothing
end

tridiag_sum(Fm1, F, Fp1, a0, am1, a, ap1) = a0 + am1*Fm1   +  a*F + ap1*Fp1


function tridiag_ldiv!( dl::AbstractVector{T}, d::AbstractVector{T}, du::AbstractVector{T}, rhs::AbstractVector{T}) where T
    n = length(rhs)
    (n == length(d) && n == length(dl) + 1 && n == length(du) + 1) || error(DimensionMismatch("All d and rhs must be of the same size and dl and du must be one element smaller"))

    @inbounds @fastmath for i = 1 : n - 1
        fact = dl[i] / d[i]           # Elimination factor
        d[i+1] -= fact * du[i]        # Update next diagonal
        rhs[i+1] -= fact * rhs[i]     # Update next RHS
    end
    
    # Back substitution
    rhs[n] = rhs[n] / d[n]
    @inbounds @fastmath  for i = n - 1 : -1 : 1
        rhs[i] = (rhs[i] - du[i] * rhs[i+1]) / d[i]
    end
    
    return rhs
end
tridiag_ldiv!(M::Tridiagonal{T,Vector{T}}, rhs::AbstractVector{T}) where T = tridiag_ldiv!(M.dl,M.d,M.du, rhs)

function column_sym_tridiag_ldiv!(rhs::AbstractVector{T}, work::AbstractVector{T}, F::AbstractVector{T}, a1::T, a0::T, a::T) where T
    n = length(rhs)
    (n == length(work) && n == length(F) ) || error(DimensionMismatch("All `F`, `rhs` and `work` must be of the same size"))
    @inbounds @fastmath for i in 1 : n 
        work[i] = a1 *  F[i] # upper diagonal indexing
        F[i] = a0 + a *  F[i]
    end  
   @inbounds @fastmath for i in 1 : n - 1                  # save original F[i]
        fact = work[i + 1] / F[i]    # subdiagonal[i] / diagonal[i]
        F[i+1] -= fact * work[i] # update next diagonal  
        rhs[i+1] -= fact * rhs[i]     # update RHS
    end
    
    rhs[n] /= F[n]
    @inbounds @fastmath for i in n-1:-1:1
        rhs[i] = (rhs[i] - work[i] * rhs[i+1]) / F[i]
    end
    
    return rhs
end