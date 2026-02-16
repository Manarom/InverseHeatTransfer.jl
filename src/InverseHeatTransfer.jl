module InverseHeatTransfer
    using LinearAlgebra, Reexport,StaticArrays
    export HeatTransferProblem
    # Write your package code here.
    include(joinpath(".","solvers", "OneDHeatTransfer.jl"))
    @reexport using .OneDHeatTransfer
    include(joinpath(".","polynomials", "PolynomialWrappers.jl"))
    @reexport using .PolynomialWrappers

    abstract type AbstractInverseProblem end
    abstract type AbstractRegularization end
    struct NoRegularization end
    abstract type AbstractCovariance end
    struct NoCovariance end    
    struct OptimizableVariable{N, DT, P,  B, V}
        p::P
        flag::B 
        lb::V
        ub::V
        is_u_bounded::Base.RefValue{Bool}
        is_l_bounded::Base.RefValue{Bool}
        function OptimizableVariable(p::ScaledPolynomial{P}; lb = nothing, ub = nothing, flag::Union{Bool, AbstractVector{Bool}, NTuple{N,Bool}} = true) where P <: AbstractPoly{N,T} where {N,T} 
                is_u_bounded = Ref(~isnothing(ub))
                is_l_bounded = Ref(~isnothing(lb))
                V = MVector{N,T}
                ub = is_u_bounded[] ? V(ub) : V(undef)
                lb = is_l_bounded[] ? V(lb) : V(undef)
                B = MVector{N,Bool}
                if isa(flag, Bool) 
                    flag_vec = B(undef)
                    fill!(flag_vec,flag)
                else
                    flag_vec = B(flag)
                end
                return new{N, T, typeof(p),  B, V}(p, flag_vec, lb, ub, is_u_bounded, is_l_bounded)
        end
    end
    parnumber(::OptimizableVariable{N}) where N = N
    coeffs(o::OptimizableVariable) = PolynomialWrappers.coeffs(o.p)
    isoptimizable(ov::OptimizableVariable) = any(ov.flag)
    refresh!(ov::OptimizableVariable, x) = any(ov.flag) ?  PolynomialWrappers.refill!(ov.p , x , ov.flag)  : nothing
    refill!(ov::OptimizableVariable, x) = PolynomialWrappers.refill!(ov.p , x)
    (ov::OptimizableVariable)(x) = ov.p(x)

    function count_lower_bound_violations(ov::OptimizableVariable{N,DT}) where {N,DT}
        ov.is_l_bounded[] || return  0
        counter = 0
        c = coeffs(ov)
        @inbounds for i in 1 : N 
           ov.flag[i] || continue
           c[i] <  ov.lb[i]  || continue     
           counter += 1
        end
        return counter
    end
    function count_upper_bound_violations(ov::OptimizableVariable{N,DT}) where {N,DT}
        ov.is_u_bounded[] || return  0
        counter = 0
        c = coeffs(ov)
        @inbounds for i in 1 : N 
           ov.flag[i] || continue
           c[i] >  ov.ub[i]  || continue     
           counter += 1
        end
        return counter

    end
    count_bound_violations(o::OptimizableVariable) = count_lower_bound_violations(o) + count_upper_bound_violations(o)
    
    function change_flag(o::OptimizableVariable; new_flag = true)
        isa(new_flag, Bool) ? fill!(o.flag, new_flag) : copyto!(o.flag, new_flag)
    end

    struct SingleInverseProblem{DT <: Number, 
                                TN , N , # TN - couples number, N - timesteps number
                                TData,
                                ProblemType <: HeatTransferProblem ,
                                CV <: AbstractCovariance, 
                                RG<: AbstractRegularization,
                                DV } 
        
        thermocouple_locations::SVector{TN , DT} # coordinates of all thermocouples
        thermocouple_indices::SVector{TN  , Int} # indices of internal thermocouples in problem TMAT  - temperature distribution matrix 
        thermocouple_values::TData # values of measured temperatures over time , number of rows  - N, 
        # number of columns must be equal to the number of locations 
        total_thickness::DT
        direct_problem::ProblemType
        covariance::CV # covariance matrix
        regularization::RG # regularization matrix
        Tdata_measured::Matrix{}
        Tdata_evaluated::DV
        residual_vector::Vector{DT} # raw residual vector
        function SingleInverseProblem(time_data ,
                                    temperatures , 
                                    initial_distribution,
                                    thermocouples_locations::AbstractVector{DT},
                                    C::OptimizableVariable,
                                    λ::OptimizableVariable, 
                                    dλdT::OptimizableVariable,
                                    thickness, 
                                    xpoints_number::Int, 
                                    tpoints_number::Int,
                                    covariance::AbstractCovariance = NoCovariance(), 
                                    regularization::AbstractRegularization= NoRegularization() 
                                    ) where DT
            
            
            issorted(time_data) || error("Time data must be sorted in ascending order")
            issorted(thermocouple_locations) || error("Thermocouple locations must be sorted in ascending order")
            NT == length(thermocouple_locations) # number of couples points
            NT < 3 && error("There should be at least three thermocouples to solve the inverse problem")
            length(time_data) == size(temperatures, 1) || error("Number of rows in temperature data should be the same as the numbe rof time points")
            
            isa(initial_distribution, VecOrMat{DT}) && length(initial_distribution) != xpoints_number && error("Number of initial distribution vector must be ")
            
            NT == size(temperatures, 2) || error("Number of thermocouple locations must 
                        be equal to the number of columns in temperatures matrix and the length of time vector ")
            
            (tmin , tmax) = extrema(time_data)
            @. time_data -= tmin
            tmax -= tmin
            
            #upper_bc_fun = Interpolations.linear_interpolation(time,)
            #=

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
            
            =#
            T_loc = SVector{NT,DT}(thermocouples_locations) # thermocouple locations 
            thickness_internal = T_loc[end] - T_loc[1]

            grid = UniformGrid(thickness_internal , tmax , Val(xpoints_number) , Val(tpoints_number))

            bc_fun_up = Interpolations.linear_interpolation(time_data , temperatures[:,1])
            bc_fun_dwn = Interpolations.linear_interpolation(time_data , temperatures[:,end])

            if isa(initial_distribution, Number)
                initT_f = InitialTFunction(Returns(initial_distribution))
            elseif isa(initial_distribution , VecOrMat)
                initT_f = InitialTFunction(linear_interpolation( collect(xrange(grid)) , initial_distribution) )
            else
                initT_f = InitialTFunction(initial_distribution)
            end

            C_f = PhysicalPropertyFunction(C)
            L_f = PhysicalPropertyFunction(λ)
            Ld_f = PhysicalPropertyFunction(dλdT)

            pr = HeatTransferProblem(C_f , L_f , Ld_f , 
                                                initT_f ,
                                                grid ,
                                                bc_fun_up ,  DirichletBC(),
                                                bc_fun_dwn , DirichletBC())
        end
    end
"""
    This type of problems include only physical properties modification, thus all problems 
has the same objects for  λ, λ' and cₚ, hence the problem can be simplified
"""
    struct PropertyInversion{TP <: Tuple, N, T}
        problems::TP
        params
        loss_vect::MVector{N,T}
    end
    function loss_function(u, ip::PropertyInversion{TP, N, T}) where {TP,N,T}
        fill_params!( u, ip)
        for (i, p_i) in enumerate(ip.problems)
            error("todo")

        end
    end
end
