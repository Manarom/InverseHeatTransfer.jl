module InverseHeatTransfer
    using LinearAlgebra , Reexport , StaticArrays , Interpolations, RecipesBase
    export HeatTransferProblem
    # Write your package code here.
    include(joinpath(".","solvers", "OneDHeatTransfer.jl"))
    @reexport using .OneDHeatTransfer
    include(joinpath(".","polynomials", "PolynomialWrappers.jl"))
    @reexport using .PolynomialWrappers

    abstract type AbstractInverseProblem end
    abstract type AbstractRegularization end
    struct NoRegularization <: AbstractRegularization end
    abstract type AbstractCovariance  end
    struct NoCovariance <: AbstractCovariance end    

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
                                ProblemType <: HeatTransferProblem ,
                                CV <: AbstractCovariance, 
                                RG <: AbstractRegularization,
                                DV } 
        
        thermocouple_locations::Vector{DT} # coordinates of all thermocouples
        thermocouple_indices::Vector{Int} # indices of internal thermocouples in problem TMAT  - temperature distribution matrix 
        thermocouple_values::Matrix{DT} 
        # values of measured temperatures over time , number of rows  - N, 
        # number of columns must be equal to the number of locations 
        total_thickness::DT
        direct_problem::ProblemType
        covariance::CV # covariance matrix
        regularization::RG # regularization matrix
        Tdata_measured::Matrix{DT}
        Tdata_evaluated::DV
        residual::Matrix{DT} # raw residual vector
        # jacobian::Matrix{}
        """
    SingleInverseProblem(
                                        time_data ::Vector{DT},
                                        temperatures ::Matrix{DT}, 
                                        initial_distribution::Union{Number, OptimizableVariable, Matrix{DT}},
                                        thermocouples_locations::AbstractVector{DT},
                                        C::OptimizableVariable,
                                        λ::OptimizableVariable, 
                                        dλdT::OptimizableVariable,
                                        thickness, 
                                        xpoints_number::Int, 
                                        tpoints_number::Int,
                                        covariance::CV = NoCovariance(), 
                                        regularization::RG = NoRegularization(),
                                        upper_flux::Union{OptimizableVariable , Nothing} = nothing,
                                        lower_flux::Union{OptimizableVariable , Nothing} = nothing,
                                        ::Type{G} = UniformGrid,
                                        thermocouple_location_relative_tolerance::Float64 = 1e-3

                                    ) where {DT , G <: AbstractGrid, CV <: AbstractCovariance, RG <: AbstractRegularization}

TBW
"""
function SingleInverseProblem(
                                        time_data ::Vector{DT},
                                        temperatures ::Matrix{DT}, 
                                        initial_distribution::Union{Number, OptimizableVariable, Matrix{DT}},
                                        thermocouples_locations::AbstractVector{DT},
                                        C::OptimizableVariable,
                                        λ::OptimizableVariable, 
                                        dλdT::OptimizableVariable,
                                        thickness, 
                                        xpoints_number::Int, 
                                        tpoints_number::Int,
                                        covariance::CV = NoCovariance(), 
                                        regularization::RG = NoRegularization(),
                                        upper_flux::Union{OptimizableVariable , Nothing} = nothing,
                                        lower_flux::Union{OptimizableVariable , Nothing} = nothing,
                                        ::Type{G} = UniformGrid,
                                        thermocouple_location_relative_tolerance::Float64 = -1.0

                                    ) where {DT , G <: AbstractGrid, CV <: AbstractCovariance, RG <: AbstractRegularization}
            
            # if fluxes are provided than the problem will be formulated with Neuman BC 
            is_upper_flux_provided = !isnothing(upper_flux) 
            is_lower_flux_provided = !isnothing(lower_flux)

            # if upper or lower heat flux is provided as an input, than we need less temperatures
            temperatures_needed = 3 - is_upper_flux_provided - is_lower_flux_provided
            issorted(time_data) || error("Time data must be sorted in ascending order")
            issorted(thermocouples_locations) || error("Thermocouple locations must be sorted in ascending order")
            NT = length(thermocouples_locations) # number of couples points
            all(Base.Fix2(<, thickness), thermocouples_locations) || error("Thermocouple locations should be smaller than the value of thickness")
            NT < temperatures_needed && error("There should be at least $(temperatures_needed) thermocouples to solve the inverse problem")
            length(time_data) == size(temperatures, 1) || error("Number of rows in temperature data should be the same as the numbe rof time points")
            
            isa(initial_distribution, VecOrMat{DT}) && length(initial_distribution) != xpoints_number && error("Number of initial distribution vector must be ")
            
            NT == size(temperatures, 2) || error("Number of thermocouple locations must 
                        be equal to the number of columns in temperatures matrix")
            

            # we need to solve the equation only in the region of interest, thus  
            # only the part of the sample is covered with grid 
            # thickness is the real thickness of the sample 
            upper_grid_coordinate = is_upper_flux_provided ? 0.0 : thermocouples_locations[1]
            lower_grid_coordinate = is_lower_flux_provided ? thickness : thermocouples_locations[end]
            thickness_internal = lower_grid_coordinate - upper_grid_coordinate
            (tmin,tmax) = extrema(time_data)
            tmin != 0.0 && (@. time_data -=tmin)
            grid = G(thickness_internal , tmax , Val(xpoints_number) , Val(tpoints_number))
            # the first and the last index of temperature columns in temperatures matrix which are used 
            first_index = is_upper_flux_provided ? 1 : 2
            last_index  = is_lower_flux_provided ? NT : NT - 1

            # here is the number of residual columns of the input data matrix which will be used for the discrepancy 
            n_residual_columns = last_index - first_index + 1
            thermocouple_indices = fill(0, (n_residual_columns,))
            rtol = thermocouple_location_relative_tolerance <= 0.0 ? 1/(2*(xpoints_number - 1)) : thermocouple_location_relative_tolerance
            
            located_inds_number = locate_indices_on_grid!(thermocouple_indices , 
                                    view(thermocouples_locations, first_index : last_index) ,
                                     grid , upper_grid_coordinate , 
                                     thickness * rtol )

            (located_inds_number != n_residual_columns) && error("Failed to attribute all thermocouple locations to the indices of grid, try to reduce the thermocouple location tolerance or the number of coordinate steps")
            
            # setting upper BC
            (bc_fun_up, bc_up_type) = if !is_upper_flux_provided 
                (Interpolations.linear_interpolation(time_data , temperatures[:,1]), DirichletBC())
            else
                (upper_flux, NeumanBC())
            end
            # setting lower BC
            (bc_fun_dwn, bc_dwn_type) = if !is_lower_flux_provided
                (Interpolations.linear_interpolation(time_data , temperatures[:,end]), DirichletBC())
            else
                (lower_flux, NeumanBC())
            end
            # setting initial temperature distribution
            x_grid = collect(xrange(grid)) 
            if isa(initial_distribution, Number) # single scalar value
                initT_f = InitialTFunction(Returns(initial_distribution), x_grid)
            elseif isa(initial_distribution , VecOrMat)
                initT_f = InitialTFunction(linear_interpolation( x_grid, initial_distribution), x_grid )
            else # provided as callable
                initT_f = InitialTFunction(initial_distribution , x_grid)
            end
            # setting physical properties
            C_f = PhysicalPropertyFunction(C)
            L_f = PhysicalPropertyFunction(λ)
            Ld_f = PhysicalPropertyFunction(dλdT)
            # setting the direct problem
            direct_problem = HeatTransferProblem(C_f , L_f , Ld_f , 
                                                initT_f ,
                                                grid ,
                                                bc_fun_up ,  bc_up_type,
                                                bc_fun_dwn , bc_dwn_type)
            t_points = tpoints(grid)
            TMAT = direct_problem.T
            Tdata_evaluated = transpose(@view TMAT[thermocouple_indices , :])# direct problem stores temperature distribution over coordinate as columns
            T_locations = collect(thermocouples_locations) # thermocouple locations  - coordinates 
            # must extract the values of measured temperatures and interpolate them of grid
            t_grid = collect(eachtime(grid))
            Tdata_measured = Matrix{DT}(undef, t_points, n_residual_columns) 
            interpolate_matrix!(Tdata_measured, time_data , temperatures , t_grid , first_index, last_index)
            residual = @. Tdata_evaluated -  Tdata_measured

            TN = temperatures_needed
            N = t_points
            ProblemType = typeof(direct_problem)
            DV = typeof(Tdata_evaluated)

            #={DT <: Number, 
                                TN , N , # TN - couples number, N - timesteps number
                                ProblemType <: HeatTransferProblem ,
                                CV <: AbstractCovariance, 
                                RG <: AbstractRegularization,
                                DV } =#
            new{DT, TN, N , ProblemType , CV , RG , DV}(        
                                                        T_locations, # thermocouple_locations - total locations including those used in BC
                                                        thermocouple_indices, # indices of thermocouples in the direct problem output matrix 
                                                        copy(temperatures), # thermocouple_values -  just copy of the input data 
                                                        thickness, # total_thickness total thickness of the sample includes the region of direct problem solution
                                                        direct_problem, #direct_problem::ProblemType direct problem solution 
                                                        covariance, # covariance matrix 
                                                        regularization, # reularization matrix 
                                                        Tdata_measured, # measured data used to evaluate the discrepancy
                                                        Tdata_evaluated, # reference to the part of temperature distribution matrix which is used for discrepancy evaluation
                                                        residual # matrix used to store the residual 
                )

        end
    end
    function fill_residual!(p::SingleInverseProblem)
        solve_problem!(p.direct_problem)
        @. p.residual = p.Tdata_evaluated - p.Tdata_measured
    end
    function interpolate_matrix!(Mout, t , M , tnew , start_col = 1, stop_col = 0)
        stop_col == 0 && (stop_col = size(M,2))
        for (i , c) in enumerate(eachcol(M)[start_col : stop_col])
            interpolator = linear_interpolation(t , c)
            c_data = @view Mout[: , i] 
            @. c_data = interpolator(tnew)
        end       
    end
    function locate_indices_on_grid!(indices_vector , locations_vector , g::AbstractGrid , zero_shift , atol )
        counter = 0
        N = length(locations_vector)
        for (i, xi) in enumerate(eachx(g))
            N <= counter && return counter
            if abs(locations_vector[counter + 1] + zero_shift - xi) <= atol
                indices_vector[counter + 1] = i
                counter += 1
            end
        end
        return counter
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
    @recipe function f(m::OptimizableVariable)
        return (m.p)
    end
end
