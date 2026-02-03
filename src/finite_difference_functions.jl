
module OneDHeatTransfer

using LinearAlgebra #,NonlinearSolvers, StaticArrays #, LaTeXStrings

using AllocCheck

const USE_FASTMATH = true
ni_err() = throw(DomainError("not implemented"))
# this types are used to dispatch on different schemes  
abstract type AbstractTimeScheme end # left-hand side part 
abstract type AbstractCoordinateScheme end # right-hand side part 
abstract type AbstractNonLinearPart end # non-linear part type

# singletones for various schemes
struct BFD1 <: AbstractTimeScheme end # first order backward difference
struct BFD2 <: AbstractTimeScheme end # second order backward difference

struct EXP <: AbstractCoordinateScheme end # fully explicit 
struct IMP <: AbstractCoordinateScheme end # fully implicit 
struct CN <: AbstractCoordinateScheme end # crank - nicolson

struct EXP_NL <: AbstractNonLinearPart end # explicit non - linear part
struct IMP_NL <: AbstractNonLinearPart end # implicit non - linear part


# boundaries conditions types 


const CentralDifference = Ref(Tridiagonal(fill(-1.0, 2),fill(0.0, 3), fill(1.0, 2))) # stores the finite difference matrix 


include("FunctionsWrappers.jl")
include("AbstractGrid.jl")
include("TridiagFunctions.jl")

# Base.iterate(g::EachTimeStep) = 

const COMMON_DOC = """

     SCHEME_NAME finite difference scheme

"""



export BFD1_EXP_EXP_EXP,BFD1_IMP_EXP_EXP,BFD1_CN_EXP_EXP

abstract type AbstractSolverScheme end
abstract type AbstractFDScheme{T,X,N,P} <: AbstractSolverScheme end
"""
Flag type FDSolverScheme{T,X,N,P} 
with parameters :
T - time derivative approximation BFD1 or BFD2
X - coordinate second derivative approximation 


"""
struct FDSolverScheme{T,X,N,P} <: AbstractFDScheme{T,X,N,P}
    FDSolverScheme(::T,::X,::N,::P) where {T <: AbstractTimeScheme,
                                         X <: AbstractCoordinateScheme,
                                         N <: AbstractNonLinearPart,
                                         P <: AbstractNonLinearPart} = new{T,X,N,P}()
end

#=
const BFD1_CN_EXP_EXP = FDSolverScheme{BFD1,CN,EXP_NL,EXP_NL}
const BFD1_IMP_EXP_EXP = FDSolverScheme{BFD1,IMP,EXP_NL,EXP_NL}
const BFD1_EXP_EXP_EXP = FDSolverScheme{BFD1,IMP,EXP_NL,EXP_NL}
=#
# generating docstrings

const SCHEME_NAMES = [:BFD1_EXP_EXP_EXP,
                      :BFD1_IMP_EXP_EXP,
                      :BFD1_CN_EXP_EXP,
                      :BFD2_IMP_EXP_EXP,
                      :BFD2_CN_EXP_EXP]
const AVAILABLE_SCHEMES = Dict{Symbol, Type{<:FDSolverScheme}}()
for d in SCHEME_NAMES
    sd = string(d)
    types_vec = split(sd,"_")
    full_name = replace(sd,"BFD1" => "first-order-backward =",
                            "EXP" => "explicit",
                            "IMP" => "implicit",
                            "CN" => "Crank-Nicolson",
                            "BFD2" =>"second-order-backward =",
                            "_" => " + ")
    @eval begin
        @doc """
            Finite difference scheme:
            time derivative scheme =  + second_order_derivative + nonlinear_part + material_properties
            $($(full_name))
        """
        const $d = FDSolverScheme{$(Symbol(types_vec[1])),
                                        $(Symbol(types_vec[2])),
                                        $(Symbol(types_vec[3]*"_NL")),
                                        $(Symbol(types_vec[4]*"_NL"))}

    end
    cur_scheme = eval(d)
    AVAILABLE_SCHEMES[d] = cur_scheme
end
function (::Type{FDSolverScheme{T,X,N,P}})() where {T <: AbstractTimeScheme,
                                         X <: AbstractCoordinateScheme,
                                         N <: AbstractNonLinearPart,
                                         P <: AbstractNonLinearPart} 
    FDSolverScheme(T(),X(),N(),P())
end

include("ExplicitSolver.jl")

allocate_tridiagonal(N::Int, T::DataType) = Tridiagonal(Vector{T}(undef, N - 1),Vector{T}(undef, N),Vector{T}(undef, N - 1))

"""
    central_finite_difference(n)

Returns central finite difference first order derivative matrix 
"""
function central_finite_difference(n)
    if size(CentralDifference[],1) != n
         CentralDifference[] = Tridiagonal(fill(-1.0, n - 1),fill(0.0, n), fill(1.0, n - 1))   
    end  
    return  CentralDifference[]
end


"""
    fill_LHS!(LHS, Fm1, F, Fp1, m, solver_scheme, problem)

Fills left-hand side of finite difference scheme inside the loop over time 

# Arguments
    - LHS - matrix to be filled 
    - Fm1 - Fm-1 lower diagonal fourier number vector (view of F[2 : M])
    - F -   Fm  - main diagonal fourier number vector (vector)
    - Fp1 - Fm+1 upper diagonal fourier number vector (view of F[1 : M - 1])
    - m - current iteration number
    - solver_scheme - see [`FDSolverScheme`](@ref)
    - problem  - PDE problem see [`HeatTransferProblem`](@ref)

"""
function fill_LHS!(LHS, Fm1, F, Fp1, m, solver_scheme, problem) 
    throw(MethodError(fill_LHS!,(LHS, Fm1, F, Fp1, m, solver_scheme, problem)))
end

"""
    fill_RHS!(b, D, Tm, Tmm1, Fm1, F, Fp1, m, phi, solver_scheme, problem)

# Arguments
    - b - righthand side vector 
    - D - finite difference matrix
    - Tm - current step temperature distribution view
    - Tmm1 - previous step temperature distribution view
    - Fm1 - Fm-1 lower diagonal fourier number vector (view of F[2 : N, m])
    - F -   Fm  - main diagonal fourier number vector (view F[:, m])
    - Fp1 - Fm+1 upper diagonal fourier number vector (view of F[1 : N - 1, m])
    - m - current iteration number
    - phi - non-linear part coefficient λ'/4λ
    - solver_scheme - see [`FDSolverScheme`](@ref)
    - problem  - PDE problem see [`HeatTransferProblem`](@ref)

"""
function fill_RHS!(b, D, Tm, Tmm1, Fm1, F, Fp1, m, phi, solver_scheme, problem)
     throw(MethodError(fill_RHS!,(b, D, Tm, Tmm1, Fm1, F, Fp1, m, phi, solver_scheme, problem)))
end

"""
    apply_bc!(dir::AbstractBCDirection, bc_type::AbstractBoundaryCondition, 
                LHS, b,  F, bc_fun ,  
                time_scheme::AbstractTimeScheme, coord_scheme::AbstractCoordinateScheme,
                t, Tm, m, dx, dt, N)

Applies boundary conditions to the finite - difference scheme matricies

# Arguments
    - LHS - lefthand side matrix reference (can be modified for some schemes)
    - b - righthand side vector, BC modifies the first of the last element accroding to the bc_fun sprcification
    - bc_fun - boundary conditions function see [`BoundaryFunction`](@ref)
    - LHS - LHS matrix, can be modified
    - F  - main diagonal fourier number vector (view F[:, m])   
    - Tm - current step temperature distribution view
    - phi - non-linear term coefficient vector
    - m -  current time iteration index
    - solver_scheme - see [`FDSolverScheme`](@ref)
    - problem  - PDE problem see [`HeatTransferProblem`](@ref)
"""
function apply_bc!(LHS, b, bc_fun,  F, Tm, phi ,   m,  solver_scheme, problem)  
     throw(MethodError(apply_bc!,(LHS, b, bc_fun,  F, Tm, phi ,   m,  solver_scheme, problem)))
end

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
thermal_conductivity(p::HeatTransferProblem{D}, T::D ) where D = p.L_f(T)
thermal_conductivity_derivative(p::HeatTransferProblem{D}, T::D ) where D = p.Ld_f(T)
heat_capacity(p::HeatTransferProblem{D}, T::D ) where D = p.C_f(T)
lower_boundary_condition(p::HeatTransferProblem{D}, t::D) where D  = p.bc_dwn(t)
upper_boundary_condition(p::HeatTransferProblem{D}, t::D) where D  = p.bc_up(t)
# BoundaryFunction{D, BF, V, <: NeumanBC , <: UpperBC}
lower_bc_type(::HeatTransferProblem{D, CF,LF,LDF,ITF, G, BCU, BCD }) where {D, CF,LF,LDF,ITF, G, BCU, BCD  <: BoundaryFunction{D, BF, V, BC_type}} where { BF, V, BC_type}  = BC_type
upper_bc_type(::HeatTransferProblem{D, CF,LF,LDF,ITF, G, BCU }) where {D, CF, LF, LDF, ITF, G, BCU <: BoundaryFunction{D, BF, V, BC_type}} where { BF, V, BC_type}  = BC_type
"""
    unified_fd_solver!( problem::HeatTransferProblem{DT, CF, LF, LDF, ITF, G, BCU, BCD, TMATtype},
                                    solver_scheme::FDSolverScheme{TS, CS, NLS, PS} = BFD1_IMP_EXP_EXP) where {DT, CF,LF,LDF,ITF, 
                                    G <:UniformGrid{N,M}, BCU, BCD, 
                                    TMATtype <: AbstractMatrix{DT},  
                                    TS <: AbstractTimeScheme,
                                    CS <: AbstractCoordinateScheme,
                                    NLS <: EXP_NL,
                                    PS <: EXP_NL} where {N, M}
                        
Unified solver for 1d heat transfer problems with various schemes

# Arguments

- problem - heat transfer problem object see [`HeatTransferProblem`](@ref)
- solver_scheme - scheme of solving see [`FDSolverScheme`](@ref) 

"""
function unified_fd_solver!( problem::HeatTransferProblem{DT, CF, LF, LDF, ITF, G, BCU, BCD, TMATtype},
                                    solver_scheme::FDSolverScheme{TS, CS, NLS, PS} = BFD1_IMP_EXP_EXP) where {DT, CF,LF,LDF,ITF, 
                                    G <:UniformGrid{N,M}, BCU, BCD, 
                                    TMATtype <: AbstractMatrix{DT},  
                                    TS <: AbstractTimeScheme,
                                    CS <: AbstractCoordinateScheme,
                                    NLS <: EXP_NL,
                                    PS <: EXP_NL} where {N, M}
        #T = Matrix{DType}(undef,N,M)# columns - distribution, rows time
        T = problem.T
        T1 = @view T[:,1]

        map!(problem.initT_f, T1, eachx(problem.grid))

        dd = problem.grid.dt/(problem.grid.dx * problem.grid.dx)#
        F, phi, LHS, D = problem.cache.F, problem.cache.phi, problem.cache.LHS, problem.cache.D
        # allocating vectors and matrices
        #F = Vector{DT}(undef,N) # properties vector
        Fm1 = @view F[2 : end] 
        Fp1 = @view F[1 : end - 1]
        # all allocations moved to the problem costructor
        #phi = Vector{DT}(undef,N) # nonlinear coefficient vector λ'/λ
        #lam = Vector{DT}(undef,N)
        #LHS = allocate_tridiagonal(N,DT) # left-hand side matrix 
        #D = central_finite_difference(N) # creates finite difference matrix for b vector evaluation
        #(C_f, L_f, Ld_f) = (,, )
        (bc_up, bc_dwn) = (problem.bc_up, problem.bc_dwn)
        #Tm = T1
        Tmm1 = T1
        # allocating left matrix
        
        is_show = N <= 10
        for m = 1 : M - 1 #% цикл по времени
            Tm = @view T[:,m] # Tm current time
            # filling current values of physical quantities
              @inbounds  for ii in 1 : N
                ti = Tm[ii]
                #abs(ti - Tmm1[ii]) >= 1e-8 || continue
                λ =  problem.L_f(ti) # λ
                #lam[ii] = λ
                F[ii] = dd * λ / problem.C_f(ti) # Fm - (dx^-2)*dt*Cp/λ
                phi[ii] = 0.25 * problem.Ld_f(ti)/λ #phi  - λ'/(4λ)
            end 

            Tmp1 = @view T[:, m + 1] # Tm+1 next time 
            !is_show || begin println("----------------------")
                            @show m
                        end
            fill_LHS!(LHS, Fm1, F, Fp1, m, solver_scheme, problem)

            fill_RHS!(Tmp1, D, Tm, Tmm1, Fm1, F, Fp1, m, phi, solver_scheme, problem)
            if is_show
                println("defore BC")
                @info LHS
                @info Tmp1
            end

            apply_bc!(LHS, Tmp1, bc_up,  F, Tm, phi ,  m,  solver_scheme, problem)

            apply_bc!(LHS, Tmp1, bc_dwn, F, Tm, phi ,  m,  solver_scheme, problem)

            if is_show
                println("after BC")
                @info LHS
                @info Tmp1
            end
            #ldiv!(Tmp1 , LHS , b) # solving LHS*T = b
            #ldiv!(LHS , b)
            tridiag_ldiv!(LHS,Tmp1)
            (TS <: BFD1) || (Tmm1 = Tm) # Tm-1 next time 
        end
   return problem
end
function evaluate_virtual_node( Fm::T, phim::T, T2::T, T0::T) where T
    return Fm * phim * (T2 - T0)^2
end

include("dirichlet_bc.jl")
include("bfd1_imp_exp_exp.jl") # fully implicit solver
include("bfd1_cn_exp_exp.jl") # crank-nicolson solver
include("bfd2_imp_exp_exp.jl") # second order backward difference
include("bfd2_cn_exp_exp.jl") # second order backward difference 
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


end