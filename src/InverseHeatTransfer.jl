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
    ov_covariance(ov::OptimizableVariable{N, DT, P, B, V} , cov , τ)  where {N, DT, P <: ScaledPolynomial, B, V}

Returns the covariance matrix of OV evaluated on vectgor τ , cov is coefficients covariance matrix 
provided externally  
"""
function  ov_covariance(ov::OptimizableVariable{N, DT, P, B, V} , cov , τ )  where {N, DT, P <: ScaledPolynomial, B, V}
            vander_matrix = ScaledPolynomials.vander(deepcopy(ov.p) , τ)
            return vander_matrix * cov * transpose(vander_matrix)
        end
        """
        create_index_mapping(optimizable)

    Used to create mapping of optimization variables vector to OptimizableVariable's Tuple
    each element of index mapping Tuple contains indices of particular optimizable variable 
    e.g. λ parameters in the input parameters vector provided by the optimizer, e.g.:

    If the optimization problem is stated for `λ` and `C` and e.g. `λ` has `3` parameters 
    and `C` has `4` parameters , than the input optimizer state vector has `7` parameters
    index mapping is `( 1:3 , 4:7)`
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

        """
        residual_length(::SingleInverseProblem{DT, TN, N} ) where {DT , TN, N}

    The length of residuals as if it is a single vector 
    """
    residual_length(::SingleInverseProblem{DT, TN, N} ) where {DT , TN, N} = N * TN
        """
        optimizable_functions_number(::SingleInverseProblem{DT, TN, N , ProblemType , CV , RG , DV, O, ON }) where {DT, TN, N , ProblemType , CV , RG , DV, O, ON }

    Total number of fucntions to be optimized 
    """
    optimizable_functions_number(::SingleInverseProblem{DT, TN, N , ProblemType , CV , RG , DV, O, ON }) where {DT, TN, N , ProblemType , CV , RG , DV, O, ON } = ON
        """
        solve_direct_problem!(p::SingleInverseProblem)

    Solves direct problem 
    """
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
    """
    optimizable_parnumber(p::SingleInverseProblem)

Total number of parameters to be modified during th eoptimization
"""
optimizable_parnumber(p::SingleInverseProblem) = sum(optimizable_parnumber, p.optimizable)
"""
    optimizable_function_names(::SingleInverseProblem{DT, TN, N , ProblemType , CV , RG , DV, O}) where {DT, TN, N , ProblemType , CV , RG , DV, O}

Returns names of functions included in the optimizable variables 
"""
optimizable_functions_names(::SingleInverseProblem{DT, TN, N , ProblemType , CV , RG , DV, O}) where {DT, TN, N , ProblemType , CV , RG , DV, O} = fieldnames(O)
    
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
        loss = 0.5 * covariance_loss(p) # applies weighted least squares (each loss is divided by the length of residual vector )
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
                    regularization = sum(p -> p.α[] * regularization_loss(p), parallel_probls.problems)
                )
    end
    function loss_distribution_matrix(parallel_probls::ParallelInverseProblems)
        return (
                    total = [evaluate_loss(p) for p in  parallel_probls.problems],
                    covariance  =  [0.5 * covariance_loss(p) for p in  parallel_probls.problems],
                    constraints =[p.ψ[] * constraints_loss(p) for p in  parallel_probls.problems],
                    regularization = [ p.α[] * regularization_loss(p) for p in parallel_probls.problems]

                )
    end
    optimizable_parnumber(p::ParallelInverseProblems) = optimizable_parnumber(first(p.problems))
    optimizable_functions_number(p::ParallelInverseProblems) = optimizable_functions_number(first(p.problems))
    optimizable_functions_names(p::ParallelInverseProblems) = optimizable_functions_names(first(p.problems))

    const ALL_REGULARIZATION_TYPES = subtypes(AbstractRegularization)
    const ALL_COVARIANCE_TYPES = subtypes(AbstractCovariance)


    """
    extract_current_solution_vector!(u , p::SingleInverseProblem)

Returns current solution state 
"""
function extract_current_solution_vector!(u , p::SingleInverseProblem) 
        copyto!(u , 
            Iterators.flatten(
                Iterators.map(fview_coeffs, p.optimizable)
            )
        )
        return u
    end
    extract_current_solution_vector!(u , p::ParallelInverseProblems) = extract_current_solution_vector!(u , first(p.problems))
    extract_current_solution_vector(p::SingleInverseProblem{DT}) where DT = extract_current_solution_vector!(Vector{DT}(undef , optimizable_parnumber(p)) , p)
    extract_current_solution_vector(p::ParallelInverseProblems) = extract_current_solution_vector(first(p.problems))
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
        σ² # estimated dispersion
        s² # sample dispersion
        sse 
        sst 
        r² # rsquared
        r²a # rsquared adjusted
        
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

            @eval struct $wrapper_type{P , T, N , M } <: AbstractStaticWrapper{P,T,N,M}
                        problem::P
                        problem_shadow::P
                        u₀::T 
                        #cache::CT
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
function discrepancy(x , p::AbstractStaticWrapper; is_specific::Bool = false)
        loss=discrepancy!(x , p.problem_shadow) 
        is_specific && (loss /= residual_length(p.problem))
        #default_state(p)
        return loss
    end

    residual(x::AbstractVector{T} , p::AbstractStaticWrapper) where {T} = residual!(
                                            Vector{T}(undef , residual_length(p.problem) ), 
                                             x , p)

    function residual!(r::AbstractVector ,   x  , p::AbstractStaticWrapper)
        discrepancy!(x , p.problem_shadow)
        copyto!(r , extract_weighted_residual_vector(p.problem_shadow ))
        #default_state(p)
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
        # default_state(p)
        return t_measured
    end

    ## Bunch of functions to evaluate the Jacobian , Hessian , sensitivity and gradient using FiniteDiff 
    # based on StaticWrapper structures, which are callable warppers around the problem restoring their 
    # state 
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
    fdif_jacobian!(J::AbstractMatrix  , srw::Union{StaticResidualWrapper , StaticEvaluatedWrapper})  = FiniteDiff.finite_difference_jacobian!(J , srw , srw.u₀  , relstep = 1e-3)

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

# this makes versions of fdif functions when calling solely on the problem , the central point of local 
# is taken from the current state of the problem 
    for f in (:fdif_gradient , :fdif_hessian , :fdif_jacobian , :fdif_sensitivity , :fdif_approximate_hessian)
        @eval $f(p::AbstractInverseProblem) = $f(p , extract_current_solution_vector(p))
        f! = Symbol(String(f)*"!")
        @eval $f!(a , p::AbstractInverseProblem) = $f!(a , extract_current_solution_vector(p) , p)
    end

# Inverse problems descriptive stats 
    """
    ip_covariance(ip::AbstractInverseProblem , u)

Function evaluates several quantities on post-processing , bu default uses approximate hessian 
J'*J but if `use_approximate_hessian` is false than uses full hessian, when there is more than 
one optimizable function this flag sticked to `true` 
"""
function ip_covariance(p::AbstractInverseProblem , u::T ; use_approximate_hessian::Bool = true ) where T <: AbstractVector

        N = residual_length(p)
        
        P = optimizable_parnumber(p)
        @assert length(u) == P "Vector  `u` length should be equal to the number of optimizable variables $(P)" # parameters number
        if use_approximate_hessian
            J = Matrix{eltype(T)}(undef, (N , P))
            fdif_jacobian!(J ,  u , p)
            return ip_approximate_covariance(J , p , u)
        else # using 

        end

    end
function ip_approximate_covariance(J::AbstractMatrix{T} , p::AbstractInverseProblem , u) where T
        OVN = optimizable_functions_number(p) # total number of optimization variables (some of them has no optimizable parameters)
        OVnames = optimizable_functions_names(p)
        N = residual_length(p)
        σ = evaluate_loss(p)/(N - optimizable_parnumber(p))
        p_i = isa(p , SingleInverseProblem) ? p : first(p.problems) # taking problem 
        (OVN > 1) && (use_approximate_hessian = true) # frocing to use approximate hessian if there is more than two optimization variables 
        out = NamedTuple{OVnames}(
            ntuple(OVN) do i
            P_i = optimizable_parnumber(p_i.optimizable[i])
            (P_i <= 0) && return nothing
            H = Matrix{T}(undef, (P_i , P_i)) # approximate jacobian 
            _J = @view J[: , p_i.index_mapper[i]] # taking only part of jacobian 
            mul!(H , transpose(_J) , _J)    
            (_regression_covariance(H , σ)... , N = N , J = _J , u = u[p_i.index_mapper[i]])
        end)
        return out
end
function ip_hessian_covariance(H , p::AbstractInverseProblem , u)
        σ² = evaluate_loss(p)
        fdif_hessian!(H  , u , p)
        N = residual_length(p)
        return (_regression_covariance(H , σ²)... , N = N ,  u=u) 
end
function _regression_covariance(H , σ²) 

    Σ = H\I 
    Cov =Σ * σ²
    s = sqrt.(diag(Cov)) 
    return (; std = s , Cov = Cov , H = H , Σ = H , σ² = σ² )
end
"""
    Type to wrap the covariance of the optimization problem 
    
    S - named tuple of statistics on different properties (names correspond to optimizables)
    ST - named tuple each name correspond to statistics type (see `ip_covariance` function)
    IsAp - true if approximate hessian was used to evaluate the covariance 
    N - number of the optimization variables 
"""
struct IPCovariance{S , ST , IsAp , N}
    stats::S
    function IPCovariance(p::AbstractInverseProblem ; use_approximate_hessian::Bool = true)
        _stats = ip_covariance(p , use_approximate_hessian = use_approximate_hessian)
        stats = filter(t->!isnothing(t) , _stats)
        new{typeof(stats) , typeof(first(stats)) , use_approximate_hessian , length(stats)}(stats)
    end
    function IPCovariance(stats::NamedTuple ; use_approximate_hessian::Bool = true)
        _stats = filter(t->!isnothing(t) , stats)
        new{typeof(_stats) , typeof(first(_stats)) , use_approximate_hessian , length(_stats)}(_stats)
    end
end
function Base.getproperty(c::IPCovariance{S , ST , IsAp , N} , p::Symbol) where {S , ST ,  IsAp , N} 
    if hasfield(IPCovariance , p)
        getfield(c , p)
    elseif hasfield(S , p)
        getfield(c.stats , p)
    elseif hasfield(ST , p)
        NamedTuple{fieldnames(S)}(
            ntuple(N) do i 
                s_i = c.stats[i]
                getfield(s_i , p)
            end
        )
    else
        error("There is no such field $(p)")
    end
end

ip_covariance(p::AbstractInverseProblem; kwargs...) = ip_covariance(p , extract_current_solution_vector(p); kwargs...)
functions_names(::IPCovariance{S}) where S = fieldnames(S)
statistics_names(::IPCovariance{S,ST}) where {S,ST} = fieldnames(ST)
has_function_name(::IPCovariance{S} , n::Symbol) where S = hasfield(S , n)

confidence_bounds(p::AbstractInverseProblem ,  τ ; kwargs...) = confidence_bounds(p , IPCovariance(p) , τ ; kwargs...)
confidence_bounds(p::ParallelInverseProblems , cov::IPCovariance ,  τ ; kwargs...) = confidence_bounds(first(p.problems) , cov , τ ; kwargs...)

"""
    confidence_bounds(p::SingleInverseProblem  , cov::IPCovariance , τ::V ) where V <: AbstractVector

Function evaluates confidence bounds for all optimizable variables associated with the inverse problem `p` 
`cov` is the inverse problem covariance object , τ is the vector of independent variable values ,  
"""
function confidence_bounds(p::SingleInverseProblem  , cov::IPCovariance , τ::V ; α::Union{Tuple , Float64} = 1.96) where V <: AbstractVector #
    # {DT , TN , N , ProblemType , CV , RG , DV , O , ON , IM}
    # {S , ST , IsAp , N}
    names = filter(n -> has_function_name(cov , n), optimizable_functions_names(p) )# names of all optimizable function 
    N = length(names)
    _α =  isa(α , Tuple) ? α : ntuple(Returns(α) , N )
    NamedTuple{names}( 
        ntuple(N) do  i
            n = names[i]
            ov = getfield( p.optimizable , n)
            c = getproperty(cov , n) # returns covariance matrix 
            confidence_bounds(ov , c.Cov , τ; α = _α[i])
        end
    )
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
        # df = (degrees_of_freedom > 0) ?  (h - degrees_of_freedom) : h # Adjust as needed with p
            p_value[i] = ccdf(Chisq(h_i), Q)
        end    

        return p_value
    end

    function sensitivity_analysis_statistics(p::AbstractInverseProblem , u) 

        J = fdif_sensitivity(p , u)
        # stats = ip_approximate_covariance(J , p , u)

        OVN = optimizable_functions_number(p) # total number of optimization variables (some of them has no optimizable parameters)
        OVnames = optimizable_functions_names(p)

        p_i = isa(p , SingleInverseProblem) ? p : first(p.problems) # taking problem 
        
        out = NamedTuple{OVnames}(
            ntuple(OVN) do i
            !is_optimizable(p_i.optimizable[i]) && return nothing 
            _J = @view J[: , p_i.index_mapper[i]] # taking only part of jacobian 
            sensitivity_analysis_statistics(_J)
        end)

        return filter(t->!isnothing(t) , out)
    end
    # 
    sensitivity_analysis_statistics(p::AbstractInverseProblem) = sensitivity_analysis_statistics(p , extract_current_solution_vector(p))
    function sensitivity_analysis_statistics(J::AbstractMatrix)
        H = transpose(J)*J
        return (
            T = t_optimality_information(H),
            D = d_optimality_information(H), 
            K = k_optimality_information(H)
        )
    end

    t_optimality_sensitivity(J) = sumsqr(J)
    t_optimality_information(H) = tr(H)
 
    d_optimality_information(H) = log(det(H))
    d_optimality_sensitivity(J) = log(det(cholesky(transpose(J)*J)))

    k_optimality_information(H) = cond(H)
    k_optimality_sensitivity(J) = cond(J)^2

    @recipe function f(m::OptimizableVariable)
        return (m.p)
    end
    include("problem_ensemble_functions.jl")
    include("hdf5_interface.jl")
    #include("tables_data.jl")
    
#########################----STATISTICS---TABLES----######################

    function autocorrelation_stats_table(p::AbstractInverseProblem; is_unweighted::Bool = false)
        s = autocorrelation_analysis(p , is_unweighted=is_unweighted)
        return _autocorrelation_stats_table(s, is_unweighted=is_unweighted)
    end
    """
        autocorrelation_stats_table(autocor_stats; is_unweighted::Bool = false)

    Creates args ready for table data creation 
    """
    function _autocorrelation_stats_table(autocor_stats; is_unweighted::Bool = false)
        N = length(autocor_stats) # number of problems 
        a1 = autocor_stats[1][1]
        fn = filter(f->isa(getfield(a1 , f) , Number) , fieldnames(typeof(a1))) # tests_names 
        M  = length(fn)
        tbl_vect = Any[]
        for (i , a_i) in enumerate(autocor_stats)
            # a_i - problem stats (vector, each element - corresponds to )
            TPnumber = length(a_i) # length of  	
            tbl_i = Matrix{Any}(undef, (TPnumber ,  M + 1) ) # i'th thermocouple table
            for tp_i in 1 : TPnumber #over couples
                a_i_t = a_i[tp_i] #
                row = Any[]
                push!(row, "P$(i) : T$(tp_i)")
                for (j , f_cur) in enumerate(fn)
                    push!(row , getfield(a_i_t, f_cur))
                end
                tbl_i[tp_i , :] = row
            end
            push!(tbl_vect , tbl_i)
        end
        col_names = Vector{String}(undef,  M + 1)
        col_names[1] = "names"
        for i in 1 : M 
            col_names[i + 1] = "$(fn[i]) "
        end	
        tbb = vcat(tbl_vect...)
        return return (tbb , ( column_labels = col_names ,  
                                title = "Autocor. tests for $(is_unweighted ? "unweighted res." : "weighted res.")"))
    end
    const AVAILABLE_STATS_TABLES =(:loss_distribution , :descriptive , :regression , :autocorrelation_weighted , :autocorrelation_unweighted ,   :sensitivity)

    function stats_table(p::AbstractInverseProblem , stats_type::Symbol; kwargs...)
        if (stats_type == :autocorrelation_weighted) || (stats_type == :autocorrelation)
            autocorrelation_stats_table(p ;is_unweighted = false, kwargs...)
        elseif stats_type == :autocorrelation_unweighted
            autocorrelation_stats_table(p ;is_unweighted = true, kwargs...)
        elseif stats_type == :descriptive
            descriptive_stats_table(p)
        elseif stats_type == :regression
            regression_stats_table(p  ; kwargs...)
        elseif stats_type == :sensitivity 
            sensitivity_stats_table(p ; kwargs...)
        elseif stats_type == :loss_distribution  
            loss_distribution_table(p ; kwargs...)
        else
            error("got unsupported stats type $stats_type ")
        end
    end

    function descriptive_stats_table(p::AbstractInverseProblem)

        ips = IPstats(p)
        fn = fieldnames(IPstats)
        N  = length(fn)
        tbl = Matrix{Any}(undef, (N , 2) )

        for i in 1:N 
            f_cur = fn[i]
            tbl[i , 1] = String(f_cur)
            tbl[i , 2] = getfield(ips , f_cur)
        end

        return (tbl ,
            (
                column_labels = ["name", "value"] ,
                title ="Simple descriptive statistics")
            )
    end

    regression_stats_table(p::AbstractInverseProblem;kwargs...) = regression_stats_table(IPCovariance(p);kwargs...)
    """
    regression_stats_table(stats , u; α = 1.96)

Regression analysis statistics table , if some of stats elements is `nothing` 
it propagates further 
"""
function regression_stats_table(stats::IPCovariance{S , ST , IsAp , N}  ; α = 1.96) where {S , ST , IsAp , N} 
    #N  = length(stats)   
    names = functions_names(stats)
    
    NamedTuple{names}(
        ntuple(N) do i 
            n = names[i]
            st = getproperty( stats , n ) 
            isnothing(st) && return nothing
            std_d = st.std
            NN = length(std_d)
            tbl = Matrix{Any}(undef, (NN , 4) )
            for ii in 1 : NN 
                tbl[ii , 1] = "$(names[i]):β$(ii)"
                tbl[ii , 2] = st.u[ii]
                tbl[ii , 3] = std_d[ii]
                tbl[ii , 4] = α * std_d[ii]
            end

            return (
                tbl ,
                (column_labels = ["name", "value" ,  "std",  "d" ] ,
                title ="Optimization variable $(names[i]) "
                )
            )
        end
    )    
    end
    """
    sensitivity_stats_table(p::AbstractInverseProblem)

Sensitivity table with optimal criteria
"""
sensitivity_stats_table(p::AbstractInverseProblem) = sensitivity_stats_table(p , extract_current_solution_vector(p))
    function sensitivity_stats_table(probs::AbstractInverseProblem , u)
        out = sensitivity_analysis_statistics(probs , u)
        tbl = Matrix{Any}(undef, (length(out) , 2))
        for (i , f_i) in enumerate(pairs(out)) 
            tbl[i , 1] = "$(first(f_i))-optimality" 
            tbl[i , 2] = last(f_i) 
        end
        
        return ( 
            tbl ,
            (column_labels = ["name"  , "value"] ,
            title ="Optimality criteria")
        )
    end
    """
    all_stats_tables(p::AbstractInverseProblem)

Returns named tuple of all avaliable stat tables data 
"""
function all_stats_tables(p::AbstractInverseProblem)
        return NamedTuple{AVAILABLE_STATS_TABLES}(
            ntuple(length(AVAILABLE_STATS_TABLES)) do i 
             stats_table(p , AVAILABLE_STATS_TABLES[i] )
        end
        )
    end
    function data_selection_table(all_data::DataSelectorsGroup)
        table_data = Matrix{Any}(undef , (length(all_data.d) , 6))
        for (i , (k , d)) in enumerate(all_data.d)
            table_data[i, 1]  = "P$(i)"
            table_data[i, 2]  = k
            table_data[i, 3] = 1e3 * DataConnector.thickness(d) # thickness is converted from m to mm 
            table_data[i, 4] = [ Pair(v,k) for (k,v) in zip(1e3 * DataConnector.sensors_locations(d) , DataConnector.selected_names(d))]
            table_data[i, 5] = DataConnector.tmin(d)
            table_data[i, 6] = DataConnector.tmax(d)
        end
        return (
                table_data , (column_labels = ["probl", "name"," h" , "Locations","tmin", "tmax"] , 
							   title ="Sample properties")
        )
    end
    function loss_distribution_table(p::AbstractInverseProblem)
        loss_table = loss_distribution_matrix(p)
        PN = isa(p , ParallelInverseProblems) ? length(p.problems) : 1
	    projects_names = ["P$(i)" for i in 1:PN]

        table_data = hcat(projects_names , [v for v in loss_table]...) 
        column_labels  =vcat("name", [String(k) for k in keys(loss_table)]...)
       return (
                table_data , 
                (title = "Loss distribution" , 
                column_labels = column_labels)
                )
    end



    ############################LCURVE_ANALYSIS############################

    
function l_curve_analysis(p::AbstractInverseProblem , alphas , solver )

	lcurve_probs = deepcopy(p)
	
	(_start, _lb, _ub)  =  fill_starting_vectors(lcurve_probs)
    r = solver(lcurve_probs)
	sols = fill(r , length(alphas))
	Threads.@sync for i in 1:length(alphas)
		 Threads.@spawn begin
			p_i = deepcopy(lcurve_probs)
			set_regularization_multiplier!(p_i , alphas[i])
            sols[i] = solver(p_i)
		end
	end
    cov_loss = similar(alphas)
	reg_loss = similar(alphas)
	total_loss = similar(alphas)
    for (i , r) in enumerate(sols) 
        discrepancy!(r.u , lcurve_probs)
        _loss_full = loss_distribution(lcurve_probs)
        cov_loss[i] = _loss_full.covariance
        reg_loss[i] = _loss_full.regularization/alphas[i]
        total_loss[i] = _loss_full.total
    end
    return (;alphas = alphas , total_loss = total_loss , cov_loss = cov_loss , reg_loss = reg_loss , sols = sols)
end
end ##end_of_module

