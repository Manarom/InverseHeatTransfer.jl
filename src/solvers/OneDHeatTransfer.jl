module OneDHeatTransfer
using LinearAlgebra
include("FunctionsWrappers.jl")
include("AbstractGrid.jl")
include(joinpath("..", "Utils","TridiagFunctions.jl"))

export HeatTransferProblem, BFD1_EXP_EXP_EXP,BFD1_IMP_EXP_EXP,BFD1_CN_EXP_EXP

struct ProblemCache{Vtype, Mtype, N, DT}
    F::Vtype
    phi::Vtype
    D::Mtype
    LHS::Mtype
    function ProblemCache(N::Int,::Type{DT}) where {DT} 
        Vtype = Vector{DT}
        F = Vtype(undef,N) # properties vector
        phi = Vtype(undef,N) # nonlinear coefficient vector λ'/λ
        LHS = allocate_tridiagonal(N,DT) # left-hand side matrix 
        Mtype = typeof(LHS)
        D = central_finite_difference(N) # creates finite difference matrix for b vector evaluation
        new{Vtype, Mtype, N, DT}(F,phi,D,LHS)
    end
end
"""
    HeatTransferProblem(C_fun, L_fun, Ld_fun, initT_fun,
                                H::D, N::Int,
                                tmax::D, M::Int,
                                bc_fun_up, upper_bc_type::AbstractBoundaryCondition,
                                bc_fun_dwn, lower_bc_type::AbstractBoundaryCondition,
                                ::Type{DT}) where {D, DT <:Number}


Formulates the problem to solve the following equation:

    (C/λ)*Tₜ= Tₓₓ + (λ'/λ)*(Tₓ)²
    T(x,0) = Tᵢ(x)

where
λ - thermal conductivity, Kg/m^3 * J/(Kg*K)    
C - thermal capacity,    C = Cp*ρ    Cp - specific heat, J/(kg*K), ρ - density, kg/m³  
Tₜ = ∂T/∂t  
Tₓ = ∂T/∂x
Tₓₓ = ∂²T/∂x² 

Dirichlet conditions:
    T(0,t) = f(t)
    T(H,t) = g(t)

Dirichlet conditions:
    Tₓ(0,t) = f(t)
    Tₓ(H,t) = g(t)

Robin conditions:
    Tₓ(0,t) = f(T)
    Tₓ(H,t) = g(T)

        HeatTransferProblem(C_fun, L_fun, Ld_fun, initT_fun,
                                H::D, N::Int,
                                tmax::D, M::Int,
                                bc_fun_up, upper_bc_type::AbstractBoundaryCondition,
                                bc_fun_dwn, lower_bc_type::AbstractBoundaryCondition,
                                ::Type{DT}) where {D, DT <:Number}
# Arguments

    - C_fun - thermal capacity , (cp*Ro) (Kg/m^3 * J/(Kg*K))
    - L_fun - thermal conductivity, W/m*K
    - Ld_fun - thermal conductivity derivative with respect to temperature
    - initT_fun - function to evaluate the initial temperature distribution  
    - H - thickness, m
    - N - number of coordinate nods
    - tmax - time interval, s
    - M - number of time nodes
    - bc_fun_up - upper boundary conditions function
    - upper_bc_type - type of upper BC: DirichletBC, NeumanBC or RobinBC
    - bc_fun_dwn - lower boundary conditions function
    - lower_bc_type - type of lower BC: DirichletBC, NeumanBC or RobinBC
    - DT - temperature data type 
"""
struct HeatTransferProblem{D, CF,LF,LDF,ITF, G, BCU, BCD,TMATtype, CacheType} # D is for temperature data type
    C_f::CF 
    L_f::LF
    Ld_f::LDF
    initT_f::ITF
    grid::G
    bc_up::BCU
    bc_dwn::BCD
    T::TMATtype
    cache::CacheType
    function HeatTransferProblem(C_fun, L_fun, Ld_fun, initT_fun,
                                H::D, N::Int,
                                tmax::D, M::Int,
                                bc_fun_up, upper_bc_type::AbstractBoundaryCondition,
                                bc_fun_dwn, lower_bc_type::AbstractBoundaryCondition,
                                ::Type{DT}) where {D, DT <:Number}

        @assert isconcretetype(DT) "Functions return type should be concrete!"

        @assert N >= 2 "Value of N must be greater than 2"
        @assert M >= 2 "Value of M must be greater than 2"
        g = UniformGrid(H,tmax,Val(N),Val(M))
        G = typeof(g)
        time_range = trange(g)


        C_f = PhysicalPropertyFunction(C_fun, nothing, DT)
        CF = typeof(C_f)
        L_f = PhysicalPropertyFunction(L_fun, nothing, DT)
        LF = typeof(L_f)
        Ld_f = PhysicalPropertyFunction(Ld_fun, nothing, DT)
        LDF = typeof(Ld_f)
        initT_f = InitialTFunction(initT_fun, time_range, DT)
        ITF = typeof(initT_f)
        TMATtype = Matrix{DT}
        T = TMATtype(undef,N,M)
        bc_up = BoundaryFunction(bc_fun_up,time_range, upper_bc_type, UPPER_BC,DT)
        BCU = typeof(bc_up)
        bc_dwn = BoundaryFunction(bc_fun_dwn, time_range,lower_bc_type, LOWER_BC, DT)
        BCD = typeof(bc_dwn)
        cache = ProblemCache(N,DT)
        CacheType = typeof(cache)
        return new{DT, CF, LF, LDF, ITF, G, BCU, BCD, TMATtype, CacheType}(C_f, L_f, Ld_f, initT_f, g, bc_up, bc_dwn, T, cache)
    end
end

tvalue(p::HeatTransferProblem, m::Int)=  tvalue(p.grid,m)
xvalue(p::HeatTransferProblem, m::Int) =  xvalue(p.grid,m)
timestep(p::HeatTransferProblem, m::Int = 1) = timestep(p.grid,m)
xstep(p::HeatTransferProblem, m::Int = 1) = xstep(p.grid,m)
xrange(p::HeatTransferProblem)= xrange(p.grid)
trange(p::HeatTransferProblem) = trange(p.grid)


thermal_diffusivity(p::HeatTransferProblem{D}, T::D ) where D = p.L_f(T)/p.C_f(T)
fourier_number(p::HeatTransferProblem{D},T, m::Int = 1) where D = thermal_diffusivity(p,T)*timestep(p,m)/(xstep(p,m)^2)
thermal_conductivity(p::HeatTransferProblem{D}, T::D ) where D = p.L_f(T)
thermal_conductivity_derivative(p::HeatTransferProblem{D}, T::D ) where D = p.Ld_f(T)
heat_capacity(p::HeatTransferProblem{D}, T::D ) where D = p.C_f(T)
lower_boundary_condition(p::HeatTransferProblem{D}, t::D) where D  = p.bc_dwn(t)
upper_boundary_condition(p::HeatTransferProblem{D}, t::D) where D  = p.bc_up(t)
# BoundaryFunction{D, BF, V, <: NeumanBC , <: UpperBC}
lower_bc_type(::HeatTransferProblem{D, CF,LF,LDF,ITF, G, BCU, BCD }) where {D, CF,LF,LDF,ITF, G, BCU, BCD  <: BoundaryFunction{D, BF, V, BC_type}} where { BF, V, BC_type}  = BC_type
upper_bc_type(::HeatTransferProblem{D, CF,LF,LDF,ITF, G, BCU }) where {D, CF, LF, LDF, ITF, G, BCU <: BoundaryFunction{D, BF, V, BC_type}} where { BF, V, BC_type}  = BC_type


"""
Bunch of functions to solve the non-linear transient heat transfer using finite difference

    a(T)Tₜ= Tₓₓ + (λ'/λ)*(Tₓ)²
    T(0,t) = f(t)
    T(H,t) = g(t)
    T(x,0) = Tᵢ(x)

    Name convention of the functions:

    lefthand side  _ second derivative _ nonlinear part _ material properties
    e.g.:
    BFD1_imp_exp_exp  => first order time derivative,
    implicit scheme for the second derivative, explicit nonlinear part, explicit physical properties

    imp - fully implicit
    exp - fully explicit
    CN - Crank-Nicolson
    BFD1 - first order backward derivative
    BFD2 - second order backward werivative

"""
OneDHeatTransfer

include("finite_difference_solver.jl")

end