

include("finite_difference_functions.jl")
using LinearAlgebra,BenchmarkTools
using Test
@testset begin 
N = 100
f = rand(N)
b = rand(N)
fm1 = @view f[2 : end]
fp1 = @view f[1 : end - 1]
M = Tridiagonal(collect(fm1), collect(f), collect(fp1))



    b0 = copy(b)
    b1 = copy(b)
    b2 = copy(b)

    mul!(b0,M,b)
    N > 5 || @show b0
    a0 = 0.0
    am1 = 1.0
    ap1 = 1.0
    a = 1.0
    OneDHeatTransfer.tridiag_mul!(b1, fm1,f, fp1, a0, am1, a, ap1)
    N > 5 ||@show b1
    OneDHeatTransfer.column_sym_tridiag_mul!(b2, f, a0, a , ap1)
    N > 5 ||@show b2

    @test norm(b1 - b0) ≈ 0 atol=1e-10
    @test norm(b2 - b0) ≈ 0 atol=1e-10

    println("benchmarking ")
    print("mul! : ")
    @btime mul!($b0,$M,$b)
    print("OneDHeatTransfer.tridiag_mul! : ")
    @btime OneDHeatTransfer.tridiag_mul!($b1,$fm1,$f, $fp1,a0, am1, a, ap1)
    print("OneDHeatTransfer.column_sym_tridiag_mul! : ")
    @btime OneDHeatTransfer.column_sym_tridiag_mul!($b,$f,a0,a,ap1) 


    y = rand(N)
    b = rand(N)
#@testset "testing muladd version " begin
     M = Tridiagonal(collect(fm1), collect(f), collect(fp1))


    q0 = copy(b)
    q1 = copy(b)
    q2 = copy(b)

    mul!(q0, M, b, 1.0, 1.0)
    N > 5 || @show q0
    a0 = 0.0
    am1 = 1.0
    ap1 = 1.0
    a = 1.0
    OneDHeatTransfer.tridiag_muladd!(q1,b,fm1,f, fp1, a0, am1, a, ap1)
    N > 5 ||@show q1


    OneDHeatTransfer.column_sym_tridiag_muladd!(q2,b, f, a0, a , ap1)
    N > 5 ||@show q2

    @test norm(q1 - q0) ≈ 0 atol=1e-10
    @test norm(q2 - q0) ≈ 0 atol=1e-10    


#end

    println("benchmarking Ab + c versions")
    print("mul!")
    @btime mul!($q0,M,$b,1.0,1.0)
    print("OneDHeatTransfer.tridiag_muladd!")
    @btime OneDHeatTransfer.tridiag_muladd!(q1,b,fm1,f, fp1, a0, am1, a, ap1)
    print("OneDHeatTransfer.column_sym_tridiag_muladd!")
    @btime OneDHeatTransfer.column_sym_tridiag_muladd!(q2,b, f, a0, a , ap1)



#@testset "solving tridiagonal tridiag_ldiv! and column_sym_tridiag_ldiv!" begin
    Fbase = rand(N)
    rhs = rand(N)
    a1 = 1.0
    a0 = 1.0
    a = 2.0
    dl = a1 * Fbase[2:end]
    du = a1 * Fbase[1:end-1]
    d = a0 .+ a*Fbase
    M = Tridiagonal(copy(dl),copy(d),copy(du))
    rhs_ldiv = copy(rhs)
    rhs_tridiag = copy(rhs)
    rhs_tridiag_from_tridiagonal = copy(rhs)
    rhs_sym_tridiag = copy(rhs)

    ldiv!(M,rhs_ldiv)
    OneDHeatTransfer.tridiag_ldiv!( copy(dl),copy(d),copy(du), rhs_tridiag)
    OneDHeatTransfer.column_sym_tridiag_ldiv!(rhs_sym_tridiag, copy(Fbase), copy(Fbase), a1, a0, a)
    OneDHeatTransfer.tridiag_ldiv!(Tridiagonal(copy(dl),copy(d),copy(du)), rhs_tridiag_from_tridiagonal)

    @test norm(rhs_ldiv .- rhs_tridiag) ≈ 0 atol=1e-10
    @test norm(rhs_ldiv .- rhs_sym_tridiag) ≈ 0 atol=1e-10
    @test norm(rhs_ldiv .- rhs_tridiag_from_tridiagonal) ≈ 0 atol=1e-10

    println("benchmarking ldiv versions")
    print("ldiv!")
    @btime ldiv!($M, $rhs_ldiv)
    print("tridiag_ldiv!")
    @btime OneDHeatTransfer.tridiag_ldiv!( $(copy(dl)),$(copy(d)),$(copy(du)), $rhs_tridiag)
    print("column_sym_tridiag_ldiv! calling on Tridiagonal")
    @btime OneDHeatTransfer.tridiag_ldiv!($(Tridiagonal(copy(dl),copy(d),copy(du))), $rhs_tridiag_from_tridiagonal)
    print("column_sym_tridiag_ldiv!")
    @btime  OneDHeatTransfer.column_sym_tridiag_ldiv!($rhs_sym_tridiag, $(copy(Fbase)), $(copy(Fbase)), a1, a0, a)

end
