




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


const COMMON_DOC = """

     SCHEME_NAME finite difference scheme

"""



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
                                    solver_scheme::FDSolverScheme{TS, CS, NLS, PS} = BFD2_IMP_EXP_EXP()) where {DT, CF,LF,LDF,ITF, 
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
        Fm1 = @view F[2 : end] 
        Fp1 = @view F[1 : end - 1]
        (bc_up, bc_dwn) = (problem.bc_up, problem.bc_dwn)
        
        Tmm1 = T1
        
        @inbounds for m = 1 : M - 1 #% цикл по времени
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

            fill_LHS!(LHS, Fm1, F, Fp1, m, solver_scheme, problem)

            fill_RHS!(Tmp1, D, Tm, Tmm1, Fm1, F, Fp1, m, phi, solver_scheme, problem)


            apply_bc!(LHS, Tmp1, bc_up,  F, Tm, phi ,  m,  solver_scheme, problem)

            apply_bc!(LHS, Tmp1, bc_dwn, F, Tm, phi ,  m,  solver_scheme, problem)

            tridiag_ldiv!(LHS,Tmp1)
            (TS <: BFD1) || (Tmm1 = Tm) # Tm-1 next time 
        end
   return problem
end
function evaluate_virtual_node( Fm::T, phim::T, T2::T, T0::T) where T
    return Fm * phim * (T2 - T0)^2
end

schemes_path = joinpath(".","schemes")
include(joinpath(schemes_path,"bfd1_exp_exp_exp.jl"))
include(joinpath(schemes_path,"dirichlet_bc.jl"))
include(joinpath(schemes_path,"bfd1_imp_exp_exp.jl")) # fully implicit solver
include(joinpath(schemes_path,"bfd1_cn_exp_exp.jl")) # crank-nicolson solver
include(joinpath(schemes_path,"bfd2_imp_exp_exp.jl")) # second order backward difference
include(joinpath(schemes_path,"bfd2_cn_exp_exp.jl")) # second order backward difference 

