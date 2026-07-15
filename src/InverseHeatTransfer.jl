module InverseHeatTransfer
    using LinearAlgebra , Reexport , StaticArrays , Interpolations, RecipesBase , Distributions
    using Unrolled
    using InteractiveUtils
    using Static
    using Accessors
    using HDF5
    using JLD2
    using Serialization
    using UUIDs 
    using StatsBase
    using FiniteDiff

    @reexport using ScaledPolynomials
    export OptimizableVariable, SingleInverseProblem

    include(joinpath(@__DIR__, "solvers", "OneDHeatTransfer.jl"))
    include(joinpath(@__DIR__, "data_utils", "DataConnector.jl"))

    @reexport using .DataConnector
    @reexport using .OneDHeatTransfer

    abstract type AbstractInverseProblem end


    abstract type AbstractCovariance  end
    struct NoCovariance <: AbstractCovariance end  
   
    const SupportedFlagType{N} = Union{Bool, AbstractVector{Bool}, NTuple{N,Bool}} where N
    const OPTIMIZABLE_VARIABLES_NAMES = (:λ , :C , :q_up , :q_dwn , :T₀ , :dλdT)

    const INVERSE_PROBLEM_HDF5_SERIALIZED_GROUPNAME = Ref("inverse_problem_serialized")
    const INVERSE_PROBLEM_HDF5_GROUPNAME = Ref("inverse_problem")
    const ALL_STATS_HDF5_GROUPNAME = Ref("statistics")

    include("optimizable_variables.jl")

    # OptimizableVariable implementatio around ScaledPolynomial
    const OVS = OptimizableVariable{N, DT, P,  B, V} where {N, DT, P <: ScaledPolynomial,  B, V}
        """
        OptimizableVariable(p::Q  ; lb::Union{Nothing , T , NTuple{N,T}, StaticVector{N,T}} = nothing, 
                                            ub::Union{Nothing , T , NTuple{N,T}, StaticVector{N,T}} = nothing, 
                                            flag::SupportedFlagType{N} = false,
                                            lb_violation_fun = < ,
                                            ub_violation_fun = > ) where {Q <: ScaledPolynomial{P}} where  P <: AbstractPoly{N,T} where {N,T}

        Special constructor for `ScaledPolynomial` approximation of the optimizable variable 
    """
    function OptimizableVariable(p::Q  ; lb::Union{Nothing , T , NTuple{N,T}, StaticVector{N,T}} = nothing, 
                                            ub::Union{Nothing , T , NTuple{N,T}, StaticVector{N,T}} = nothing, 
                                            flag::SupportedFlagType{N} = false,
                                            lb_violation_fun = < ,
                                            ub_violation_fun = > ) where {Q <: ScaledPolynomial{P}} where  P <: AbstractPoly{N,T} where {N,T} 
                
                is_u_bounded = Ref(~isnothing(ub))
                is_l_bounded = Ref(~isnothing(lb))

                is_u_single_number = isa(lb , Number)
                is_l_single_number = isa(ub , Number)

                V = MVector{N,T}

                _ub = (!is_u_bounded[] || is_u_single_number) ? P(V(undef)) : P(V(ub))  
                _lb = (!is_l_bounded[] || is_l_single_number) ? P(V(undef)) : P(V(lb)) 

                is_u_single_number && fill!(_ub , ub)
                is_l_single_number && fill!(_lb , lb)

                if (is_u_bounded[] && is_l_bounded[]) 
                    count_violations(ScaledPolynomials.coeffs(_ub),
                    ScaledPolynomials.coeffs(_lb), lb_violation_fun) > 0 && error("Lower boundary $(_lb) is higher than $(_ub)") 
                end    
                B = MVector{N,Bool}
                if isa(flag, Bool) 
                    flag_vec = B(undef)
                    fill!(flag_vec,flag)
                else
                    flag_vec = B(flag)
                end
                return OptimizableVariable(T, p, 
                                        flag_vec, _lb, _ub,
                                        is_u_bounded, is_l_bounded,
                                        lb_violation_fun,
                                        ub_violation_fun)
        end
        # implementing `OptimizableVariable` interface
        coeffs(o::OVS) = ScaledPolynomials.coeffs(o.p)

        lb_coeffs(ov::OVS)  =  ScaledPolynomials.coeffs(ov.lb)
        ub_coeffs(ov::OVS)  =  ScaledPolynomials.coeffs(ov.ub)

        derivative!(ov_der::OVS, ov::OVS) = ScaledPolynomials.derivative!(ov_der.p, ov.p)
        """
        create_index_mapping(optimizable)

    Used to create mapping of optimization variables vector to OptimizableVariable's Tuple
    """
    function create_index_mapping(optimizable)
                cursor = 1 # updated variables counter 
                v = []
                for ov in optimizable
                    if !is_optimizable(ov) 
                        push!(v , cursor : cursor  - 1)
                        continue
                    end
                    n = optimizable_parnumber(ov)
                    push!(v , cursor : cursor + n - 1)
                    cursor += n
                end    
                return Tuple(v)
        end
    # const POSSIBLE_TAGS = (:lam, :C, )
        """
            Type to store the single inverse problem (simplest case of the problem )
            parameters of the type:

                DT  - type of temperature data 
                TN - number of thermocouples involved in residual  (number of columns of redual matrix)
                N - number of time points 
                ProblemType - direct problem type (see `HeatTransferProblem` for details)
                CV - covarinace type 
                RG - regularization type 
                DV - evalauted data (view of full results matrix )
                O - a Tuple of optimizable variables (to change the optimizabel variables the problem should be recreated)
                ON - total number of optimizable variables 
                IM - Tuple of index mapping with the same number of elements as the optimizable, conatins indices ranges which must be used 
                to fill the parameters of the optimizable variables from the extrenal vector (see `update_all_optimizables!`)
        """
        struct SingleInverseProblem{DT <: Number, 
                                    TN , N , # TN - couples number, N - timesteps number
                                    ProblemType <: HeatTransferProblem ,
                                    CV <: AbstractCovariance, 
                                    RG <: AbstractRegularization,
                                    DV, 
                                    O <: NamedTuple, # optimizable variable iterator 
                                    ON, # number of optimizable variables (variable which can possibly be optimized) 
                                    IM <: Tuple # index mapping for optimizable variables parameters in optimization variables single vector 
                                    } <: AbstractInverseProblem
            
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
            optimizable::O # this field stores the iterable object over all variables to be optimizaed
            index_mapper::IM # stored indices ranges for optimizable variables parameters vector 
            α::Base.RefValue{DT} # regularization multiplier
            ψ::Base.RefValue{DT}  # constraints violation multiplier if contraint violation is added to the loss function 
            include_constraints_violation_to_loss::Base.RefValue{Bool} # if this flag is true, conatraints violation are added to the loss function
@doc raw"""     
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

            Units: time - seconds, coordinate - meters, temperature - Celsius

        # Input arguments:

        time_data ::Vector{DT}, - vector of times in seconds
        temperatures ::Matrix{DT},  - measured temperatures in oC
        initial_distribution::Union{OptimizableVariable, Number , VecOrMat{DT}}, - starting temperature distribution in oC, 
        can be provided in several ways:
            - as an optimizable variable (in this case it will be optimized)
            - as a fixed number (in this case initial distribution assumed constant)
            - 
        thermocouples_locations::AbstractVector{DT}, - thermocouple locations in m
        C, - volumetric heat capacity (callable object ) , C(T) returns Cₚ⋅ρ where ρ - density in kg/m³ and Cₚ is specific heat J/(kg⋅ᵒC)
        λ, - thermal conductivity (callable object ),  λ(T) - returns thermal conductivity in W/m*K
        dλdT, - thermal conductivity derivative 
        thickness::Number, - total thickness of the sample in m 
        xpoints_number::Int = 200, - number of coordinate points 
        time_points_number::Union{Int,Nothing} = 2000;  - number of time points 
        covariance::CV = NoCovariance(),  - covariance type, see [`ALL_COVARIANCE_TYPES`](@ref)
        regularization::RG = NoRegularization(), - regularization type, see [`ALL_REGULARIZATION_TYPES`](@ref) 
        upper_flux::Union{OptimizableVariable , Nothing} = nothing, 
        lower_flux::Union{OptimizableVariable , Nothing} = nothing,
        grid_type::Type{G} = UniformGrid ,
        thermocouple_location_relative_tolerance::Float64 = -1.0 ,
        alpha::Number = 1e-3,
        psi::Number = 1e-3 ,
        include_constraints_violation_to_loss::Bool =false

    Setting upper_flux and lower_flux  to nothing  automatically enforces Dirichlet BC, in this case temperatures 
    data on thermocouples with lowest and highest coordinate are taken as BC and the sample thickness of the dicrect 
    problem assumed to be the difference between this coordinates. 
    """
    function SingleInverseProblem(
                                            time_data ::Vector{DT},
                                            temperatures ::Matrix{DT}, 
                                            initial_distribution::Union{OptimizableVariable, Number , VecOrMat{DT}},
                                            thermocouples_locations::AbstractVector{DT},
                                            C,
                                            λ, 
                                            dλdT,
                                            thickness::Number, 
                                            xpoints_number::Int = 200, 
                                            time_points_number::Union{Int,Nothing} = 2000;
                                            covariance::CV = NoCovariance(), 
                                            regularization::RG = NoRegularization(),
                                            upper_flux::Union{OptimizableVariable , Nothing} = nothing,
                                            lower_flux::Union{OptimizableVariable , Nothing} = nothing,
                                            grid_type::Type{G} = UniformGrid ,
                                            thermocouple_location_relative_tolerance::Float64 = -1.0 ,
                                            alpha::Number = 1e-3,
                                            psi::Number = 1e-3 ,
                                            include_constraints_violation_to_loss::Bool =false
                                        ) where {DT , G <: AbstractGrid, CV <: AbstractCovariance, RG <: AbstractRegularization}
                
                # if fluxes are provided than the problem will be formulated with Neuman BC 
                is_upper_flux_provided = !isnothing(upper_flux) 
                is_lower_flux_provided = !isnothing(lower_flux)

                # if upper or lower heat flux is provided as an input, than we need less temperatures
                temperatures_needed = 3 - is_upper_flux_provided - is_lower_flux_provided # number of sensorces need to solve the inverse problem

                issorted(time_data) || error("Time data must be sorted in ascending order")
                issorted(thermocouples_locations) || error("Thermocouple locations must be sorted in ascending order")
                NT = length(thermocouples_locations) # number of couples points including those used in BC formulation 
                all(Base.Fix2(<=, thickness), thermocouples_locations) || error("Thermocouple locations should be smaller than the value of thickness")
                (NT < temperatures_needed) && error("There should be at least $(temperatures_needed) thermocouples to solve the inverse problem")
                (length(time_data) == size(temperatures, 1)) || error("Number of rows in temperature data should be the same as the numbe rof time points")
                
                isa(initial_distribution, VecOrMat{DT}) && length(initial_distribution) != xpoints_number && error("Number of initial distribution vector must be ")
                
                NT == size(temperatures, 2) || error("Number of thermocouple locations must 
                            be equal to the number of columns in temperatures matrix")
                

                if isa(dλdT, OptimizableVariable)
                    change_flag!(dλdT , new_flag = false )
                else
                    dλdT = OptimizableVariable(dλdT)
                end         
                # we need to solve the equation only in the region of interest, thus  
                # only the part of the sample is covered with grid 
                # thickness is the real thickness of the sample 
                upper_grid_coordinate = is_upper_flux_provided ? 0.0 : thermocouples_locations[1]
                lower_grid_coordinate = is_lower_flux_provided ? thickness : thermocouples_locations[end]
                thickness_internal = lower_grid_coordinate - upper_grid_coordinate
                (tmin , tmax) = extrema(time_data)
                @. time_data -= tmin # shifting time data to make it starting from zero
                tmax = tmax - tmin
                tpoints_number = isnothing(time_points_number) ? length(time_data) : time_points_number
                grid = G(thickness_internal , tmax , Val(xpoints_number) , Val(tpoints_number))
                # the first and the last index of temperature columns in temperatures matrix which are used 
                first_index = is_upper_flux_provided ? 1 : 2
                last_index  = is_lower_flux_provided ? NT : NT - 1

                # here is the number of residual columns of the input data matrix which will be used for the discrepancy 
                n_residual_columns = last_index - first_index + 1
                thermocouple_indices = fill(0, (n_residual_columns,))

                rtol = thermocouple_location_relative_tolerance <= 0.0 ? 1/(2*(xpoints_number - 1)) : thermocouple_location_relative_tolerance
                
                located_inds_number = locate_indices_on_grid!(thermocouple_indices , 
                                                                thermocouples_locations[first_index : last_index],
                                                                grid , 
                                                                upper_grid_coordinate , 
                                                                thickness * rtol )

                (located_inds_number != n_residual_columns) && error("Failed to attribute all thermocouple locations to the indices of grid, try to reduce the thermocouple location tolerance or the number of coordinate steps")
                
                # setting upper BC
                (bc_fun_up, bc_up_type) = if !is_upper_flux_provided 
                    (Interpolations.linear_interpolation(time_data , temperatures[:,1] , extrapolation_bc=Line()), DirichletBC())
                else
                    (upper_flux, NeumanBC())
                end
                # setting lower BC
                (bc_fun_dwn, bc_dwn_type) = if !is_lower_flux_provided
                    (Interpolations.linear_interpolation(time_data , temperatures[:,end] , extrapolation_bc=Line()), DirichletBC())
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

                #TN = temperatures_needed
                TN = n_residual_columns
                N = t_points
                ProblemType = typeof(direct_problem)
                DV = typeof(Tdata_evaluated)

                # all possibly optimizable variables are arranged into named tuple 

                #optimizable =(;λ = λ , C = C, dλdT = dλdT) # q_up = upper_flux , q_dwn = lower_flux , T₀ = initial_distribution, )
                
                optimizable =(; (k => v for (k, v) in zip(OPTIMIZABLE_VARIABLES_NAMES, (λ, C, upper_flux, lower_flux, initial_distribution, dλdT)) if isa(v, OptimizableVariable))...)
                index_mapper = create_index_mapping(optimizable)
                O = typeof(optimizable)
                ON = length(optimizable)
                IM = typeof(index_mapper)

                obj = new{DT, TN, N , ProblemType , CV , RG , DV, O, ON , IM}(        
                                                            T_locations, # thermocouple_locations - total locations including those used in BC
                                                            thermocouple_indices, # indices of thermocouples in the direct problem output matrix 
                                                            copy(temperatures), # thermocouple_values -  just copy of the input data 
                                                            thickness, # total_thickness total thickness of the sample includes the region of direct problem solution
                                                            direct_problem, #direct_problem::ProblemType direct problem solution 
                                                            covariance, # covariance matrix 
                                                            regularization, # reularization matrix 
                                                            Tdata_measured, # measured data used to evaluate the discrepancy
                                                            Tdata_evaluated, # reference to the part of temperature distribution matrix which is used for discrepancy evaluation
                                                            residual, # matrix used to store the residual
                                                            optimizable,
                                                            index_mapper,
                                                            Ref(alpha),
                                                            Ref(psi),
                                                            Ref(include_constraints_violation_to_loss) 
                    )
                fill_covariance_cache!(obj)
                return obj
            end
        end

function SingleInverseProblem(
                                data_selector :: DataConnector.DataSelector,
                                C,
                                λ, 
                                dλdT,
                                xpoints_number::Int = 200, 
                                time_points_number::Union{Int,Nothing} = 2000;kwargs...)
                (   
                    time_data,
				    temperatures, 
				    initial_distribution,
				    thermocouples_locations,
				    thickness
                ) = DataConnector.combine_selected_data(data_selector)       
                
                inds = sortperm(thermocouples_locations)
                
        return SingleInverseProblem(time_data , temperatures[:,inds] , 
                                        initial_distribution, thermocouples_locations[inds] ,
                                        C,
                                        λ, 
                                        dλdT, 
                                        thickness; kwargs...)
        end
        """
        fill_residual!(p::SingleInverseProblem)

    Function  solves the direct problem and refills the resiaduals matrix
    """
    function fill_residual!(p::SingleInverseProblem)
            @. p.residual = p.Tdata_evaluated - p.Tdata_measured
            return nothing
    end

    residual_length(::SingleInverseProblem{DT, TN, N} ) where {DT , TN, N} = N * TN

    solve_direct_problem!(p::SingleInverseProblem) = solve_problem!(p.direct_problem)


        """
        update_all_optimizables!(p::SingleInverseProblem, x_vector::AbstractVector)

    Function updates all optimizable variables with respect to the flags vectors
    """
    function update_all_optimizables!(p::SingleInverseProblem{DT, TN, N , ProblemType , CV , RG , DV, O, ON , IM}, 
        x::AbstractVector) where {DT, TN, N , ProblemType , CV , RG , DV, O, ON , IM}

            ntuple(Val(ON)) do i 
                modify!(p.optimizable[i] , x , p.index_mapper[i])
            end
            
            is_λ_optimizable(p) && modify_λ_derivative!(p)
            return nothing
        end
#= 
looks like version using ntuple is faster than `Unrolled`
    @unroll function _modify_optimizables(optimizables , index_mapper::D{V} , x) where V
        @unroll for i  in 1 : V
                (ov, r) = (optimizables[i] , index_mapper[i])
                #modify!(ov, view(x , r))
                modify!(ov, x , r)
            end
    end
=#
    optimizable_parnumber(p::SingleInverseProblem) = sum(optimizable_parnumber, p.optimizable)

    is_λ_optimizable(p::SingleInverseProblem) = haskey(p.optimizable,:λ)

    function modify_λ_derivative!(p::SingleInverseProblem) 
        derivative!(p.optimizable.dλdT, p.optimizable.λ)
    end

    """
    fill_starting_vector(p::SingleInverseProblem{DT}) where DT

Scans all optimizable variables and returns the tuple with `(;x₀ - starting vector,
lb - lower boundaries vector, ub - upper boundaries vector)``  for the optimization

"""
function fill_starting_vectors(p::SingleInverseProblem{DT}) where DT
        v = Vector{DT}()
        lb = Vector{DT}()
        ub = Vector{DT}()
        cursor = 1 # updated variables counter 
        for ov in p.optimizable

            !is_optimizable(ov) && continue
            n = optimizable_parnumber(ov)
            cur_length = length(v) + n

            resize!(v  , cur_length )
            resize!(lb , cur_length )
            resize!(ub , cur_length )

            _v  = view(v , cursor : cursor + n - 1)
            _lb = view(lb, cursor : cursor + n - 1)
            _ub = view(ub, cursor : cursor + n - 1)

            ilb = is_lower_bounded(ov)
            iub = is_upper_bounded(ov)

            _d , _l , _u = fview_coeffs(ov), fview_lb_coeffs(ov) , fview_ub_coeffs(ov)

            if  ilb && iub
               @. _v = 0.5 * (_u + _l)
               copyto!(_lb , _l)
               copyto!(_ub , _u)
            elseif ilb
                copyto!(_v , _l)
                copyto!(_lb , _l)
                fill!(_ub , DT(Inf))
            elseif iub
                copyto!(_v , _u)
                copyto!(_ub , _u)
                fill!(_lb , DT(-Inf))
            else
                copyto!(_v , _d)
                fill!( _ub , DT( Inf) )
                fill!( _lb , DT(-Inf) )
            end

            cursor += n
        end
        return (; x₀ = v , lb = lb , ub = ub)
    end

    function extract_current_solution_vector!(u , p::SingleInverseProblem{DT , TN, N, ProblemType , CV , RG , DV, O , ON}) where {DT , TN, N, ProblemType , CV , RG , DV, O , ON}
        copyto!(u , Iterators.flatten(

        )
        )
    end
    """
    constraints_violation_loss(p::SingleInverseProblem)

Evaluates constraints violation part of discrepancy
"""
constraints_loss(p::SingleInverseProblem) = sum(constraints_loss , p.optimizable)

 """
    regularization_loss(::SingleInverseProblem{DT,TN,N,P,CV,RG}) where {DT , TN , N , P,CV,RG <: NoRegularization}

No regularization
"""
regularization_loss(::SingleInverseProblem{DT,TN,N,P,CV,RG}) where {DT , TN , N , P,CV,RG <: NoRegularization} = zero(DT)
    """
    regularize(::FiniteDifferenceRegularization , p::SingleInverseProblem{DT})

Regularization with a finite difference matrix returns ``Σᵢ [xᵀDᵀDx /(<Δx>² nᵢ)]ᵢ`` the 
summation is over all `OptimizableVariable` in problem, here x is a parameters vector
for the `OptimizableVariable`, ``<Δx> = (max(x) - min(x))^2``

When using together with bernstein polynomials this regularization reduces the `steepness` of the output
function

"""
regularization_loss( p::SingleInverseProblem{DT , TN , N , P , CV , RG ,  DV ,  O , ON}) where {DT , TN , N , P,CV, RG <: FiniteDifferenceRegularization,  DV ,  O, ON} = sum(finite_difference_regularization_loss , p.optimizable)

regularization_loss( p::SingleInverseProblem{DT ,TN , N , P , CV , RG , DV , O , ON}) where {DT , TN , N , P,CV, RG <: FixedDiagonalRegularization,  DV ,  O, ON} = sum(fixed_diagonal_regularization_loss , p.optimizable)

include("covariances.jl")


    """
    discrepancy(x , p::SingleInverseProblem{DT}) where DT

Evaluates the weighted least-sqaure discrepancy of the corresponding inverse problem 
fills parameters -> solves heat transfer problem -> updates residuals -> evaluates total loss
"""
function discrepancy!(x , p::SingleInverseProblem{DT}) where DT
        update_all_optimizables!(p , x) # refreshes the values of parameters without solving the direct problem 
        solve_direct_problem!(p) # solves the direct problem 
        fill_residual!(p) # fills residual matrix 
        return evaluate_loss(p)
    end

function set_regularization_multiplier!(s::SingleInverseProblem{DT} , α::DT) where DT 
        s.α[] = α
        return nothing
    end
    """
    evaluate_loss(p :: SingleInverseProblem{DT}) where DT

Function evaluates scalar discrepancy for the current set of parameters 
"""
function evaluate_loss(p :: SingleInverseProblem{DT}) where DT
        loss = covariance_loss(p) # applies weighted least squares (each loss is divided by the length of residual vector )
        p.include_constraints_violation_to_loss[] && (loss += p.ψ[] * constraints_loss(p)) # adds constraints loss to the main discrepancy (if they are needed)
        loss += p.α[] * regularization_loss(p)
        return loss
    end
    function interpolate_matrix!(Mout, t , M , tnew , start_col::Int = 1, stop_col::Int = 0)
        
        stop_col <= 0 && (stop_col = size(M,2))
        iter_step = start_col <= stop_col ? 1 : -1
        for (i , c) in enumerate(eachcol(M)[start_col : iter_step : stop_col])
            interpolator = linear_interpolation(t , c , extrapolation_bc=Line())
            c_data = @view Mout[: , i] 
            @. c_data = interpolator(tnew)
        end       
    end

    function locate_indices_on_grid!(indices_vector , locations_vector , g::AbstractGrid , zero_shift , atol )
        counter = 0
        N = length(locations_vector)
        for (i, xi) in enumerate(eachx(g))
            N <= counter && return counter
            if abs(locations_vector[counter + 1] - zero_shift - xi) <= atol
                indices_vector[counter + 1] = i
                counter += 1
            end
        end
        return counter
    end
"""
    This type of problems include only physical properties modification, thus all problems 
has the same objects for  λ, λ' and cₚ, hence the problem can be simplified.
"""
    struct ParallelInverseProblems{TP <: Tuple, N, T} <: AbstractInverseProblem
        problems::TP
        function ParallelInverseProblems(probls::SingleInverseProblem{DT}...) where DT
            N = length(probls)
            new{typeof(probls) , N , DT}(probls )
        end
    end

    function residual_length(pp::ParallelInverseProblems{TP , N}) where {TP , N}
        return sum( ntuple(N) do i
            residual_length(pp.problems[i]) 
        end  
        )
    end
    fill_starting_vectors(pp::ParallelInverseProblems) = fill_starting_vectors(pp.problems[1])

    """
    discrepancy!(x , pp::ParallelInverseProblems{TP, N, T}) where {TP , N , T}

"""
function discrepancy!(x , pp::ParallelInverseProblems{TP, N}) where {TP , N}
        return sum(
            ntuple( N ) do i 
                discrepancy!(x , pp.problems[i])
            end
            )# /N  total discrepancy divided by the problems number 
    end
 function evaluate_loss(pp::ParallelInverseProblems{TP, N}) where {TP , N }
     return sum(evaluate_loss , pp.problems)/N # dividing by the number of problems 
 end
    set_regularization_multiplier!(pp::ParallelInverseProblems , val) = foreach(pp.problems) do p 
                                    set_regularization_multiplier!(p , val)
                                end

    function loss_distribution(p::SingleInverseProblem )
        return (
                    total = evaluate_loss(p),
                    covariance  = covariance_loss(p),
                    constraints = p.ψ[] * constraints_loss(p),
                    regularization = p.α[] * regularization_loss(p)
                )
    end
    function loss_distribution(parallel_probls::ParallelInverseProblems)
        return (
                    total = sum(evaluate_loss , parallel_probls.problems),
                    covariance  = sum(covariance_loss,  parallel_probls.problems),
                    constraints =sum(constraints_loss,  parallel_probls.problems),
                    regularization = sum(p->p.α[] * regularization_loss(p), parallel_probls.problems)
                )
    end
    function loss_distribution_matrix(parallel_probls::ParallelInverseProblems)
        return (
                    total = [evaluate_loss(p) for p in  parallel_probls.problems],
                    covariance  = [covariance_loss(p) for p in  parallel_probls.problems],
                    constraints =[p.ψ[] * constraints_loss(p) for p in  parallel_probls.problems],
                    regularization = [ p.α[] * regularization_loss(p) for p in parallel_probls.problems]

                )
    end
    optimizable_parnumber(p::ParallelInverseProblems) = optimizable_parnumber(first(p.problems))
    const ALL_REGULARIZATION_TYPES = subtypes(AbstractRegularization)
    const ALL_COVARIANCE_TYPES = subtypes(AbstractCovariance)


        """
        extract_residual_vector(p::SingleInverseProblem)

    retutrns reference to the residual vector 
    """
    #extract_residual_vector(p::SingleInverseProblem) = Iterators.flatten(p.residual)
    #extract_residual_vector(p::ParallelInverseProblems{D,N}) where {D,N}

    for (field_name , func_name) in zip((:residual, :Tdata_measured , :Tdata_evaluated), (:extract_residual_vector , :extract_measured_vector , :extract_evaluated_vector))
        str = String(field_name)
        @eval $func_name(p::SingleInverseProblem) = Iterators.flatten(getfield(p , Symbol($str)))
        @eval function $func_name(p::ParallelInverseProblems{D,N}) where {D,N}
            return Iterators.flatten(
                ntuple(N) do i 
                    getfield(p.problems[i] , Symbol($str))
                end
            )
        end

    end


    function extract_weighted_residual_vector(p::SingleInverseProblem)
        return ResidualIterator(Val(true) , p)
    end
    function extract_weighted_residual_vector(p::ParallelInverseProblems{D,N}) where {D,N}
            return Iterators.flatten(
                ntuple(N) do i 
                    ResidualIterator(Val(true) , p.problems[i])
                end
            )
    end

    struct IPstats{N ,P}
        σ2 # estimated dispersion
        s2 # sample dispersion
        sse 
        sst 
        r2 # rsquared
        r2a # rsquared adjusted
        
        function IPstats(y , r , N::Int, P::Int=1)
            sse = sumsqr(r)
            σ2  = sse/(N - P)
            m = StatsBase.mean(y)
            sst = sumsqr(y , m)
            s2 = sst/(N - 1)
            r2 = 1 - sse/sst 
            r2a = 1 - (sse/(N - 1)) * (N - P)/sst
            new{N , P}( σ2 , s2 , sse , sst , r2 , r2a )
        end
    end
    function IPstats(p::AbstractInverseProblem)
        IPstats(extract_measured_vector(p) , extract_weighted_residual_vector(p) , residual_length(p) , optimizable_parnumber(p) )
    end
    function sumsqr(itr, μ::T =0.0) where T 
            s = zero(T)
            for t in itr
                s+=(t - μ)^2
            end
            return s
    end

    """
        Wrappers for finite difference methods of derivative evaluation. Each of AbstractStaticWrapper automatically
    returns to its initial state after optimization variables change
    P - inverse problem type 
    T - type of optimization parameter vector (input vector)
    N - total number of residual points 
    M - length of the optimization variables vector 
    """
   abstract type AbstractStaticWrapper{P , T , N , M} end


    for wrapper_type in (:StaticDiscrepancyWrapper , :StaticResidualWrapper , :StaticEvaluatedWrapper)
        @eval struct $wrapper_type{P , T, N , M} <: AbstractStaticWrapper{P,T,N,M}
            problem::P
            problem_shadow::P
            u₀::T 
        end       
    end
    function (::Type{ASW})(p::P , u₀::T ) where ASW <: AbstractStaticWrapper where {P <: AbstractInverseProblem , T <: AbstractVector } 
        M = length(u₀)
        N = residual_length(p)
        ASW{P , T , N , M}(p , deepcopy(p) , u₀)
    end
    """
    (p::StaticDiscrepancyWrapper)(x)

Callable obj for discrepancy value , after evaluation returns problem to its previous state  
"""
(p::StaticDiscrepancyWrapper)(x)  = discrepancy(x , p)
    """
    (p::StaticResidualWrapper)(x)

Vector function for inverse problem residuals 
"""
(p::StaticResidualWrapper)(x) = residual(x , p  )
    """
    (p::StaticResidualWrapper)(r , x)

Can be used in-place to fill residauls vector 
"""
(p::StaticResidualWrapper)(r , x) = residual!(r , x , p  )

    """
    (p::StaticEvaluatedWrapper)(x)

Returns evaluated temperature distribution, can be used for sensitivity analysis 
"""
(p::StaticEvaluatedWrapper)(x) = evaluated(x , p  )
(p::StaticEvaluatedWrapper)(r , x) = evaluated!(r , x , p  )

    """
    default_state(p::StaticProblemWrapper)

ReturnsStaticDiscrepancyWrapper to its initial state 
"""
function default_state(p::AbstractStaticWrapper)
        # update_all_optimizables!(p.problem_shadow , p.u₀) # refreshes the values of parameters without solving the direct problem 
        # solve_direct_problem!(p.problem_shadow) # solves the direct problem 
        # fill_residual!(p.problem_shadow ) # fills residual matrix 
        discrepancy!(p.u₀ , p.problem_shadow)
    end
    """
    discrepancy(p::StaticProblemWrapper , x)

Evaluates scalar discrepancy function on input vector 'x' but 
if `is_specific` is true (default) returns total discrepancy divided 
by the total number of residual points 
"""
function discrepancy(x , p::AbstractStaticWrapper; is_specific::Bool = true)
        loss=discrepancy!(x , p.problem_shadow) 
        !is_specific && (loss *= residual_length(p.problem))
        default_state(p)
        return loss
    end

    residual(x::AbstractVector{T} , p::AbstractStaticWrapper) where {T} = residual!(
                                            Vector{T}(undef , residual_length(p.problem) ), 
                                             x , p)

    function residual!(r::AbstractVector ,   x  , p::AbstractStaticWrapper)
        discrepancy!(x , p.problem_shadow)
        copyto!(r , extract_weighted_residual_vector(p.problem_shadow ))
        default_state(p)
        return r
    end


    function evaluated(x::AbstractVector{T} , p::AbstractStaticWrapper) where {T}
        return evaluated!(
                        Vector{T}(undef , residual_length(p.problem) ), 
                        x , p)
    end
    function evaluated!(t_measured::AbstractVector ,   x  , p::AbstractStaticWrapper)
        discrepancy!(x , p.problem_shadow)
        copyto!(t_measured , extract_evaluated_vector(p.problem_shadow ))
        default_state(p)
        return t_measured
    end

    """
    fdif_hessian( p::AbstractInverseProblem , u::T) where T

Hessian matrix using FiniteDiff package
"""
function fdif_hessian( p::AbstractInverseProblem , u::T) where T 
        M = optimizable_parnumber(p)
        @assert length(u)==M "Length of u must be the same as the number of th optimizable parameters $(M)"
        H = Matrix{eltype(T)}(undef, (M , M))
        fdif_hessian!(H ,  u , p)
        return H
    end
    fdif_hessian!(H::AbstractMatrix  , u , p::AbstractInverseProblem) = fdif_hessian!(H , StaticDiscrepancyWrapper(p , u))
    fdif_hessian!(H , spw::StaticDiscrepancyWrapper ) = FiniteDiff.finite_difference_hessian!(H  , spw , spw.u₀ )

    """
    fdif_gradient( p::AbstractInverseProblem , u::T) where T

Evaluates the gradient using `FiniteDiff` package
"""
function fdif_gradient( p::AbstractInverseProblem , u::T) where T 
        M = optimizable_parnumber(p)
        @assert length(u)==M "Length of u must be the same as the number of th optimizable parameters $(M)"
        g = Vector{eltype(T)}(undef, M)
        fdif_gradient!(g ,  u , p)
        return g
    end
    """
    fdif_gradient!(g::AbstractVector , p::AbstractInverseProblem , u)

In-place version of finite difference gradient evaluation 
"""
fdif_gradient!(g::AbstractVector  , u , p::AbstractInverseProblem) = fdif_gradient!(g , StaticDiscrepancyWrapper(p , u))
fdif_gradient!(g , spw::StaticDiscrepancyWrapper ) = FiniteDiff.finite_difference_gradient!(g  , spw , spw.u₀ )
        """
        fdif_jacobian( p::AbstractInverseProblem , u::T) where T

    Jacobian matrix of the problem using FiniteDiff package 
    """
    function fdif_jacobian( p::AbstractInverseProblem , u::T) where T 
            M = optimizable_parnumber(p)
            N = residual_length(p)
            @assert length(u)==M "Length of u must be the same as the number of th optimizable parameters $(M)"
            J = Matrix{eltype(T)}(undef, (N , M))
            fdif_jacobian!(J , u , p)
            return J
        end
    fdif_jacobian!(J::AbstractMatrix , u , p::AbstractInverseProblem ) = fdif_jacobian!(J , StaticResidualWrapper(p , u))
    fdif_jacobian!(J::AbstractMatrix  , srw::Union{StaticResidualWrapper , StaticEvaluatedWrapper})  = FiniteDiff.finite_difference_jacobian!(J , srw , srw.u₀)

    """
    fdif_sensitivity( p::AbstractInverseProblem , u::T) where T

Function to evaluate the Jacobian of evaluated temperature distributions for sensitivity analysis (∇Ţcalculated)
Standard fdif_jacobian evaluates ∇r (r is the weighted residual vector )
"""
function fdif_sensitivity( p::AbstractInverseProblem , u::T) where T 
            M = optimizable_parnumber(p)
            N = residual_length(p)
            @assert length(u)==M "Length of u must be the same as the number of th optimizable parameters $(M)"
            J = Matrix{eltype(T)}(undef, (N , M))
            fdif_sensitivity!(J ,  u , p)
            return J
        end
    fdif_sensitivity!(J::AbstractMatrix  , u , p::AbstractInverseProblem) = fdif_jacobian!(J , StaticEvaluatedWrapper(p , u))


    # approximate hessian functions 
    """
    fdif_approximate_hessian!(H, p::AbstractInverseProblem , u::T) where T

Function for approximate Hessian J'J , J is weighted residuals Jacobian 
"""
function fdif_approximate_hessian!(H , u::T ,  p::AbstractInverseProblem) where T 
        J = fdif_jacobian(p , u)
        mul!(H , transpose(J) , J)
        return H
    end
    fdif_approximate_hessian( p::AbstractInverseProblem , u::T) where T = fdif_approximate_hessian!(Matrix{eltype(T)}(undef , ntuple(_->length(u) , 2)) , u , p) 
# Inverse problems descriptive stats 
    """
    ip_covariance(ip::AbstractInverseProblem , u)

Function evaluates several quantities on post-processing
"""
function ip_covariance(p::AbstractInverseProblem , u::T; use_approximate_hessian::Bool = true ) where T <: AbstractVector

        N = residual_length(p)
        σ = evaluate_loss(p)
        P = optimizable_parnumber(p)
        @assert length(u) == P "Vector  `u` length should be equal to the number of optimizable variables $(P)" # parameters number
        J = Matrix{eltype(T)}(undef, (N , P))
        fdif_jacobian!(J ,  u , p)
        H = Matrix{eltype(T)}(undef, (P , P)) # approximate jacobian 
        if use_approximate_hessian
            mul!(H , transpose(J) , J)
        else
            fdif_hessian!(H  , u , p)
        end
        #=try 
            cholesky!(H)
        catch  er 
            warning(er)
        end=#
        Σ = H\I 
        Cov =Σ * σ
        s = @view Cov[diagind(Cov)]

        return (;std = s , Cov = Cov , H = H , Σ =H , σ = σ , N = N , J = J)
    end

    """
    autocorrelation_analysis(p::ParallelInverseProblems{TP , N}; is_unweighted::Bool = false) where {TP , N}


Evaluates descriptive statistics on weighted or unweighted residuals , if `is_unweighted` unweighted 
residual vector is used (returns Tuple of Vector of Named tuples )
"""
autocorrelation_analysis(p::ParallelInverseProblems{TP , N}; is_unweighted::Bool = false) where {TP , N}= ntuple(N) do i 
                                                                                    autocorrelation_analysis(p.problems[i] , is_unweighted = is_unweighted)
                                                                                end
    function autocorrelation_analysis(p::SingleInverseProblem{DT , TN, N} ; is_unweighted::Bool = false) where {DT, TN, N}
        residuals_iterator = is_unweighted ?  eachcol(p.residual) : eachcol(extract_weighted_residual_vector(p))
        r, cr =ntuple(_-> Vector{DT}(undef, N) , 2)
        lgs = collect(0:(N - 1)) # autocorrelation lags 

        stat_eltype =  NamedTuple{(:dubin_watson , :integral_test , :autocor), Tuple{DT, DT , Vector{DT}}} 

        stats = Vector{stat_eltype}(undef , TN)

        for (i , w) in enumerate(residuals_iterator)
            copyto!(r  , w)
            autocor!(cr , r ,  lgs)
            stats[i] =  (
                            dubin_watson = dubin_watson(r) ,  
                            integral_test = integral_cor_test(cr)  , 
                            autocor = copy(cr) 
                        )
        end
        return stats
    end
    function dubin_watson(y)
        ssqr = sumsqr(y)
        s = zero(eltype(y))
        for i in 2:length(y)
            s+=(y[i] - y[i-1])^2
        end
        return s/ssqr
    end
    function integral_cor_test(cor)
        return sumsqr(view(cor , 2 : length(cor)))/sumsqr(cor)
    end
    function ljungbox(r::AbstractVector{D}) where D
        ""
        n = length(r)
        h = 2:Int(round(log(n)))
        acor = autocor(r , 1:maximum(h))
        
        p_value = Vector{D}(undef,length(h))
        
        for (i, h_i) in enumerate(h)
            Q = zero(D)
            for k = 1 : h_i
                Q +=  (acor[k]^2)/(n - k)
            end
            Q *= (n - 1) * (n + 1)
        #df = (degrees_of_freedom > 0) ?  (h - degrees_of_freedom) : h # Adjust as needed with p
            p_value[i] = ccdf(Chisq(h_i), Q)
        end    

        return p_value
    end

    sensitivity_analysis_statistics(p::AbstractInverseProblem , u) = sensitivity_analysis_statistics(fdif_sensitivity(p , u) )
    # 
    function sensitivity_analysis_statistics(J)
        H = transpose(J)*J
        return (
            T = t_optimality_information(H),
            D = d_optimality_information(H), 
            K = k_optimality_information(H)

            )
    end

    t_optimality_sensitivity(J) = sumsqr(J)
    t_optimality_information(H) = sum(diag(H))
 
    d_optimality_information(H) = log(det(H))
    d_optimality_sensitivity(J) = log(det(cholesky(transpose(J)*J)))


    k_optimality_information(H::AbstractMatrix) = cond(H)
    k_optimality_sensitivity(J::AbstractMatrix) = cond(J)^2


    @recipe function f(m::OptimizableVariable)
        return (m.p)
    end
    include("problem_ensemble_functions.jl")
    include("hdf5_interface.jl")
    
end ##end_of_module

