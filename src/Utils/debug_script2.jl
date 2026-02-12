
using LinearAlgebra
include("TridiagFunctions.jl")

N = 3
F = [11.0,22,33]
rhs = [4.0,5,12]
a1 = 1.0
a0 = 1.0
a = 2.0

dl = a1 * F[2:end]
du = a1 * F[1:end-1]
d = a .+ a0*F

M = Tridiagonal(copy(dl),copy(d),copy(du))
rhs_ldiv = copy(rhs)
rhs_tridiag = copy(rhs)
rhs_sym_tridiag = copy(rhs)

#@show M
ldiv!(M,rhs_ldiv)
#@show M
tridiag_ldiv!( copy(dl),copy(d),copy(du), rhs_tridiag)
a1
a0
a
thomas_d!(rhs_sym_tridiag, copy(F), copy(F), a1, a0, a)

norm(rhs_ldiv .- rhs_tridiag)
norm(rhs_ldiv .- rhs_sym_tridiag)