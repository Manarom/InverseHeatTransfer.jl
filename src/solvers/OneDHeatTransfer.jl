module OneDHeatTransfer
        using LinearAlgebra, Interpolations

        include("function_wrappers.jl")
        include("abstract_grid.jl")
        include( "TridiagFunctions.jl")
        
        export HeatTransferProblem , BFD1_EXP_EXP_EXP , BFD1_IMP_EXP_EXP , BFD1_CN_EXP_EXP , solve_problem!
        export DirichletBC, NeumanBC , temperature_field
        export PhysicalPropertyFunction, BoundaryFunction, InitialTFunction
        export AbstractGrid, AbstractGridIterator, UniformGrid, eachx, eachtime, trange, xrange, tpoints, xpoints , timestep

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

    """
    HeatTransferProblem(C_f::CF, 
                        L_f::LF, 
                        Ld_f::LDF, 
                        initT_f::ITF,
                        grid::G,
                        bc_fun_up, upper_bc_type::AbstractBoundaryCondition,
                        bc_fun_dwn, lower_bc_type::AbstractBoundaryCondition) where {CF<:PhysicalPropertyFunction{DT},
                                            LF<:PhysicalPropertyFunction{DT},
                                            G <: AbstractGrid{N,M,D},
                                            LDF<: PhysicalPropertyFunction{DT},
                                            ITF <: InitialTFunction{DT}} where {D <: Number, DT <:Number, N, M}

Internal constructor accepts type-wrappers for input functions, ensures type stability
"""
function HeatTransferProblem(C_f::CF, 
                                        L_f::LF, 
                                        Ld_f::LDF, 
                                        initT_f::ITF,
                                        grid::G,
                                        bc_fun_up, upper_bc_type::AbstractBoundaryCondition,
                                        bc_fun_dwn, lower_bc_type::AbstractBoundaryCondition) where {CF<:PhysicalPropertyFunction{DT},
                                                           LF<:PhysicalPropertyFunction{DT},
                                                           G <: AbstractGrid{N,M,D},
                                                           LDF<: PhysicalPropertyFunction{DT},
                                                           ITF <: InitialTFunction{DT}} where {D <: Number, DT <:Number, N, M} 
                                                           # DT - type returned by physical properties functions 
                                                           # D - type of grid 

                    time_range = trange(grid)
                    TMATtype = Matrix{DT}
                    T = TMATtype(undef,N,M)
                    bc_up = BoundaryFunction(bc_fun_up,time_range, upper_bc_type, UPPER_BC,DT)
                    BCU = typeof(bc_up)
                    bc_dwn = BoundaryFunction(bc_fun_dwn, time_range,lower_bc_type, LOWER_BC, DT)
                    BCD = typeof(bc_dwn)
                    cache = ProblemCache(N,DT)
                    CacheType = typeof(cache)

                    return new{DT, CF, LF, LDF, ITF, G, BCU, BCD, TMATtype, CacheType}(C_f, L_f, Ld_f, initT_f, grid, bc_up, bc_dwn, T, cache)
            end
        end
    """
    HeatTransferProblem(C_fun, L_fun, Ld_fun, initT_fun,
                                    H::D, N::Int,
                                    tmax::D, M::Int,
                                    bc_fun_up, bc_fun_dwn, ::Type{DT},
                                    upper_bc_type::AbstractBoundaryCondition=DirichletBC(),
                                    lower_bc_type::AbstractBoundaryCondition=DirichletBC()) where {D, DT <:Number}

        
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
            - bc_fun_dwn - lower boundary conditions function
            - DT - temperature data type (default Float64)
            - upper_bc_type - type of upper BC: DirichletBC, NeumanBC or RobinBC (default DirichletBC)
            - lower_bc_type - type of lower BC: DirichletBC, NeumanBC or RobinBC (default DirichletBC)

"""
function HeatTransferProblem(C_fun, L_fun, Ld_fun, 
                                    initT_fun,
                                    H::D, N::Int,
                                    tmax::D, M::Int,
                                    bc_fun_up, bc_fun_dwn, 
                                    ::Type{DT} = Float64,
                                    upper_bc_type::AbstractBoundaryCondition=DirichletBC(),
                                    lower_bc_type::AbstractBoundaryCondition=DirichletBC(),
                                    ::Type{GridType} = UniformGrid) where {D, DT <:Number , GridType <: AbstractGrid}
                
                @assert isconcretetype(DT) "Functions return type should be concrete!"
                @assert N >= 2 "Value of N must be greater than 2"
                @assert M >= 2 "Value of M must be greater than 2"
                g = UniformGrid(H,tmax,Val(N),Val(M))
                time_range = trange(g)
                C_f = PhysicalPropertyFunction(C_fun, nothing, DT)
                L_f = PhysicalPropertyFunction(L_fun, nothing, DT)
                Ld_f = PhysicalPropertyFunction(Ld_fun, nothing, DT)
                initT_f = InitialTFunction(initT_fun, time_range, DT)
                return HeatTransferProblem(C_f , L_f , Ld_f ,  initT_f ,
                                            g ,  
                                            bc_fun_up, upper_bc_type , 
                                            bc_fun_dwn, lower_bc_type)
            end
        """
    copy_physics(p::HeatTransferProblem{DT, CF, LF, LDF, ITF, G},                                     
                            initT_fun,
                            H::D, tmax::D, 
                            bc_fun_up, bc_fun_dwn ) where {D, DT, CF, LF, LDF, ITF, G}

This function creates HeatTransferProblem with the same physical properties, but new geometry, time, init distribution and BC's
"""
function copy_physics(p::HeatTransferProblem{DT, CF, LF, LDF, ITF, G},                                     
                            initT_fun,
                            H::D, tmax::D, 
                            bc_fun_up, bc_fun_dwn ) where {D, DT, CF, LF, LDF, ITF, G}
            
            grid = G(H,tmax)
            time_range = trange(grid)
            initT_f = InitialTFunction(initT_fun, time_range, DT)
            BCU = upper_bc_type(p)
            BCD = lower_bc_type(p)
            return HeatTransferProblem(p.C_f, p.L_f, p.Ld_f,
                                        initT_f, grid,
                                        bc_fun_up, BCU() , 
                                        bc_fun_dwn, BCD())
            

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
    temperature_field(p::HeatTransferProblem)

Returns temperature interpolator: T(x,t) , here x is coordinate, t  - time 
"""
function temperature_field(p::HeatTransferProblem)
            #return interpolate((collect(xrange(p)), collect(trange(p))), p.T , Gridded(Linear()))
            return scale(interpolate(p.T,    BSpline(Linear())), xrange(p.grid), trange(p.grid))
        end


all_grid_interpolators(p::HeatTransferProblem) =(T = temperature_field(p), 
                                                 Tₓ = temperature_gradient(p),
                                                 Tₓₓ = temperature_laplacian(p))
"""
    temperature_gradient(p::HeatTransferProblem)

returns the interpolator of ∂T/∂x(x,t) on grid of the problem `p`
uses central finite difference to evaluate the in-grid points and 
second order backward differences at the boundaries
"""
function temperature_gradient(p::HeatTransferProblem)
    return scale(interpolate(first_order_x_derivative(p) ,    BSpline(Linear())), xrange(p.grid), trange(p.grid))
    #return interpolate((collect(xrange(p)), collect(trange(p))), first_order_x_derivative(p), Gridded(Linear()))
end
"""
    temperature_laplacian(p::HeatTransferProblem)

returns the interpolator of ∂²T/∂x²(x,t) on grid of the problem `p`
"""
function temperature_laplacian(p::HeatTransferProblem)
    return scale(interpolate(first_order_x_derivative(p) ,    BSpline(Linear())), xrange(p.grid), trange(p.grid))
end
first_order_x_derivative(p::HeatTransferProblem) = first_order_x_derivative!(similar(p.T),p)
        """
    first_order_x_derivative!(Tx::TMATtype , p::HeatTransferProblem{DT, CF, LF,
                                         LDF, ITF, G, BCU, BCD, TMATtype})  where {DT, CF, 
                                            LF, LDF, ITF, 
                                            G <: UniformGrid, BCU, BCD, TMATtype}

Function evaluates ∂T/∂x from te solution and returns it as a matrix
"""
function first_order_x_derivative!(Tx::TMATtype , p::HeatTransferProblem{DT, CF, LF,
                                         LDF, ITF, G, BCU, BCD, TMATtype})  where {DT, CF, 
                                            LF, LDF, ITF, 
                                            G <: UniformGrid, BCU, BCD, TMATtype}
    grid = p.grid                           
    Nx = OneDHeatTransfer.xpoints(grid)
    T_mat = p.T
    # first derivatives internal points
    dx = OneDHeatTransfer.xstep(grid, 1)
    @inbounds for i in 2 : Nx - 1
        
        Tp1 = @view T_mat[i + 1 , :]
        Tm1 = @view T_mat[i - 1 , :]
        Txc = @view Tx[i , :]
        @. Txc = (Tp1 - Tm1)/ (2dx)
    end
    # second order backward finite difference
    T1, T2, T3 = @view(T_mat[1 , :]), @view(T_mat[2 , :]), @view(T_mat[3 , :])
    @. Tx[1, :] = (-3 * T1 + 4 * T2 - T3) / (2 * dx)
    
    TN, TNm1, TNm2 = @view(T_mat[Nx , :]), @view(T_mat[Nx - 1 , :]), @view(T_mat[Nx - 2 , :])
    @. Tx[Nx, :] = (3 * TN - 4 * TNm1 + TNm2) / (2 * dx)
    
    return Tx
end
"""
    second_order_x_derivative(p::HeatTransferProblem)

Returns the second oreder derivative matrix of the same size as the temperature matrix `p.T`
on grid.
"""
second_order_x_derivative(p::HeatTransferProblem) = second_order_x_derivative!(similar(p.T),p)
"""
    second_order_x_derivative!(Txx::TMATtype , p::HeatTransferProblem{DT, CF, LF,
                                         LDF, ITF, G, BCU, BCD, TMATtype})  where {DT, CF, 
                                            LF, LDF, ITF, 
                                            G <: UniformGrid, BCU, BCD, TMATtype}

Fills the matrix of second derivatives of the same size as temperature distribution matrix inside the 
`HeatTransferProblem` object applying the central derivative

!!! In current implementation boundary value are just the same as the adjacent
"""
function second_order_x_derivative!(Txx::TMATtype , p::HeatTransferProblem{DT, CF, LF,
                                         LDF, ITF, G, BCU, BCD, TMATtype})  where {DT, CF, 
                                            LF, LDF, ITF, 
                                            G <: UniformGrid, BCU, BCD, TMATtype} 

    grid = p.grid                           
    Nx = OneDHeatTransfer.xpoints(grid)
    T_mat = p.T
    # first derivatives internal points
    dx = OneDHeatTransfer.xstep(grid, 1)
    @inbounds for i in 2 : Nx - 1
        Tp1 = @view T_mat[i + 1 , :]
        Tm1 = @view T_mat[i - 1 , :]
        Ti = @view T_mat[i  , :]
        Txxc = @view Txx[i , :]
        @. Txxc = (Tp1 -2*Ti + Tm1)/ (dx^2)
    end
    
    Txx[1, :] .= Txx[2, :]
    Txx[Nx, :] .= Txx[Nx - 1, :]

    return Txx
end
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
            CN -  Crank-Nicolson
            BFD1 - first order backward derivative
            BFD2 - second order backward derivative

        """
        OneDHeatTransfer
        include("finite_difference_solver.jl")
        solve_problem!(p::HeatTransferProblem) = unified_fd_solver!(p)
end