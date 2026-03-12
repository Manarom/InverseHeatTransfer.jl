module InverseHeatTransfer
    using LinearAlgebra , Reexport , StaticArrays , Interpolations, RecipesBase
    using Unrolled
    using InteractiveUtils
    import FunctionWrappers
    
    export OptimizableVariable, SingleInverseProblem
    # Write your package code here.
    include(joinpath(".","solvers", "OneDHeatTransfer.jl"))
    @reexport using .OneDHeatTransfer
    include(joinpath(".","polynomials", "PolynomialWrappers.jl"))
    @reexport using .PolynomialWrappers

    abstract type AbstractInverseProblem end
    abstract type AbstractRegularization end

    struct NoRegularization <: AbstractRegularization end
    struct  FiniteDifferenceRegularization <: AbstractRegularization end


    abstract type AbstractCovariance  end
    struct NoCovariance <: AbstractCovariance end  
   
 

    const SupportedFlagType{N} = Union{Bool, AbstractVector{Bool}, NTuple{N,Bool}} where N
    #=struct BinaryPredicate{D}  
        f::FunctionWrappers.FunctionWrapper{Bool,Tuple{D,D}} 
    end
    
    (f::BinaryPredicate{D})(x::D,y::D) where D = f.f(x,y)=# #this was slow

"""
    wrappes any modifiable and callcable variable of type P, which has parameters accessabel by `coeffs` function
    and can be bounded by `lb` and `ub` constraints, constraints violation can be checked  with  `lb_violation_fun`
    and `ub_violation_fun` functions. The parameters of `ub` and `lb` objects should be accessable by `lb_coeffs` and 
    `ub_coeffs`, there also should be a function to evaluate the derivative of `lb` and `ub` with respect to their 
    parameters.
    
        Currently implemented for `ScaledPolynomial` from `PolynomialWrappers`
        The following interface should be implemented to make it work for any particular type
    coeffs(::OV{N,DT,P}) where {N,DT,P} 
    lb_coeffs(::OV{N,DT,P})  
    ub_coeffs(::OV{N,DT,P}) where {N,DT,P}  
    # not necessary
    derivative!(::OV, ::OV)
"""
    struct OptimizableVariable{N, DT, P,  B, V, FL, FU}
        p::P # any type which must be callable
        flag::B 
        lb::V
        ub::V
        is_u_bounded::Base.RefValue{Bool} # flag if upper bounded 
        is_l_bounded::Base.RefValue{Bool}
        lb_violation_fun::FL
        ub_violation_fun::FU
        lb_violation::MVector{N,DT}
        ub_violation::MVector{N,DT}
        """
    OptimizableVariable(::Type{DT}, p::P, flag::B, lb::V, ub::V, 
                                 iu::Ref{Bool}, il::Ref{Bool}, 
                                 f1::F1, f2::F2) where {P, B <: MVector{N,Bool}, V, F1, F2} where {N,DT}

    Wrapper interface for some callable object of type `P `which can be mutated by index `flag`

    This structure can be bounded from the `top` and from the `bottom` using some objects `ub` and `lb` , both of the same type `V`
Variables `iu` and `il` are just flags which can be used to turn `on` and `off` this boundaries, boundaries are compared 
to the current state of the varible `p` using `f1` and `f2` functions, which must take two argument of `P` and `V` types 
and return a single value.

"""
    function OptimizableVariable(::Type{DT}, p::P, flag::B, 
                                 lb::V, ub::V, 
                                 iu::Ref{Bool}, il::Ref{Bool}, 
                                 f1::F1, 
                                 f2::F2) where {P, B <: MVector{N,Bool}, V, F1, F2} where {N,DT}
            # Здесь можно добавить проверки размеров, если нужно
            lb_violation = MVector{N,DT}(fill(zero(DT) ,  N))
            ub_violation = MVector{N,DT}(fill(zero(DT) ,  N))
            new{N, DT, P, B, V, F1, F2}(p, flag, lb, ub, iu, il, f1, f2, lb_violation, ub_violation)
        end
    end
     # default methods
    const OV = OptimizableVariable

    (ov::OV)(x) = ov.p(x)


is_optimizable(_) = false
total_parnumber(::OV{N}) where N = N
is_optimizable(ov::OV) = any(ov.flag)
is_lower_bounded(ov::OV) = ov.is_l_bounded[]
is_upper_bounded(ov::OV) = ov.is_u_bounded[]
optimizable_parnumber(ov::OV) = sum(ov.flag)
change_flag!(o::OV; new_flag = true) = isa(new_flag, Bool) ? fill!(o.flag, new_flag) : copyto!(o.flag, new_flag)
    """
    refill!(ov::OV, x)

Refills all coefficients from another vector 
"""
refill!(ov::OV, x) = begin 
        copyto!(coeffs(ov), x)
        return nothing
    end
    """
    modify!(ov::OV, x)

Modifies coefficients which are marked as adjustable by ov.flag
If x is empty or ov.flag has no true values does nothing
"""
function modify!(ov::OV{N , DT}, x) where {N , DT}
        (isempty(x) || !is_optimizable(ov)) && return nothing
        optimizable_parnumber(ov) != length(x) && error("incorrect x size") 
        #_x = coeffs(ov)
        #_f = ov.flag
        counter = 1
        @inbounds for i = 1 : N 
            if ov.flag[i]
                coeffs(ov)[i] = x[counter]
                counter += 1
            end    
        end
        # copyto!(fview_coeffs(ov), x)
        return nothing
    end
"""
    modify!(ov::OV{N , DT}, x , r) where {N , DT}

Function fills parameters of `ov` marked by ov.flag from 
x[r], r can be a vector of indices or indices range
"""
function modify!(ov::OV{N , DT}, x , r) where {N , DT}
        n = optimizable_parnumber(ov)
        nr = length(r)
        nr == 0 && return nothing
        (n != nr || nr > length(x)) && error("incorrect x size") 
        counter = 1
        @inbounds for i = 1 : N 
            if ov.flag[i]
                coeffs(ov)[i] = x[r[counter]]
                counter += 1
            end    
        end
        # copyto!(fview_coeffs(ov), x)
        return nothing
    end
    count_violations(a, b, f) = count(f(i,k) for (i,k) in zip(a,b))
    count_lower_bound_violations(ov::OV) = (is_lower_bounded(ov) && is_optimizable(ov)) ? count_violations(fview_coeffs(ov),fview_lb_coeffs(ov), ov.lb_violation_fun) : 0
    count_upper_bound_violations(ov::OV)= (is_upper_bounded(ov) && is_optimizable(ov)) ? count_violations(fview_coeffs(ov),fview_ub_coeffs(ov), ov.ub_violation_fun) : 0
    count_bound_violations(o::OV) = count_lower_bound_violations(o) + count_upper_bound_violations(o)

    """
    constraints_loss(ov::OV{N , DT}) where {N,DT}

Evaluates scalar loss due to the `OptimizableVariable` constraints violation 
The value of loss is proportional to the square of the difference between the constraint
value and the actual value of coefficients, normalized to the span of the box
If the optimization variable is box constraint, thus having bot `lb` and `ub` vectors
limiting the possible range of coefficients , constraints loss are evaluated as 

`Σᵢ[(xᵢ - lbᵢ)/spanᵢ]² + Σⱼ[(xⱼ - ubⱼ)/spanⱼ]²` where `x` is the optimizable 
variabel parameters `spanᵢ = ubᵢ - lbᵢ` is the box width, and `i`and `j` are 
the indices of coordinates  which are marked as optimizable and violate 
lower or upper constraints respectively

"""
function constraints_loss(ov::OV{N , DT}) where {N,DT}
        !is_optimizable(ov) && return zero(DT)
        ilb = is_lower_bounded(ov)
        iub = is_upper_bounded(ov)
        !iub && !ilb && return zero(DT)
        s_i = zero(DT)
        @inbounds for i in 1 : N
                if ov.flag[i]
                    x_i = coeffs(ov)[i]
                    ub_i, lb_i = ub_coeffs(ov)[i], lb_coeffs(ov)[i]
                    span = if ilb && iub
                        span = abs(ub_i -  lb_i)
                    elseif ilb
                        lb_i
                    else
                        ub_i    
                    end    
                    if ilb && ov.lb_violation_fun(x_i , lb_i)
                        ov.lb_violation[i] = ((x_i -  lb_i)/span)^2.0
                        s_i += ov.lb_violation[i] 
                    else
                        ov.lb_violation[i] = zero(DT)
                    end
                    if iub && ov.ub_violation_fun(x_i , ub_i)
                        ov.ub_violation[i] = ((x_i -  ub_i)/span)^2.0
                        s_i += ov.ub_violation[i]
                    else
                        ov.ub_violation[i] = zero(DT)
                    end
                end
            end
        return s_i/N
    end
    """
    finite_difference_regularization_loss(ov::OptimizableVariable{N,DT}) where {N,DT}

Evaluates loss addition due to Tikhonov's regularization  `xᵀDᵀDx/N` with regularizing matrix `D`
is a finite difference matrix here `x` is `ALL` coefficients vector (not only those which are 
supposed to be modified by flag) 

When using together with Bernstein polynomial forces function to be more monotonical
"""
function finite_difference_regularization_loss(ov::OV{N,DT}) where {N,DT}

        !is_optimizable(ov) && return zero(DT)
        _x = coeffs(ov)
        s_i = zero(DT)

        (lv, hv) = (_x[1], _x[1])
        
        @inbounds @simd for i in 1 : N - 1
            x_ip = _x[i  + 1]
            hv = max(hv , x_ip)
            lv = min(lv , x_ip)
            Δ = x_ip - _x[i]
            s_i += Δ * Δ
        end
        diff = hv - lv
        s_i *= (abs(diff) > DT(1e-16)) ? DT(0.25) / (N * (diff^2)) : zero(DT)
        return s_i   
    end

    fview_coeffs(ov::OV) = view(coeffs(ov), ov.flag)
    fview_lb_coeffs(ov::OV) = view(lb_coeffs(ov), ov.flag)
    fview_ub_coeffs(ov::OV) = view(ub_coeffs(ov), ov.flag)
    fview_lb_violation(ov::OV)  = view(ov.lb_violation, ov.flag)
    fview_ub_violation(ov::OV)  = view(ov.ub_violation, ov.flag)

    extract_all_params(ov::OV) = copy(coeffs(ov))
    extract_optimizable_params(ov::OV) = copy(fview_coeffs(ov))
    # interface
    # necessary
    coeffs(::OV{N,DT,P}) where {N,DT,P}  = error("OptimizableVariable wrappers is not implemenented for $(P) type")
    lb_coeffs(::OV{N,DT,P})  where {N,DT,P}  = error("OptimizableVariable wrappers around $(P) is not implemenented")
    ub_coeffs(::OV{N,DT,P}) where {N,DT,P}  = error("OptimizableVariable wrappers around $(P) is not implemenented")
    # not necessary
    derivative!(::OV, ::OV)  = error("OptimizableVariable wrappers around is not implemented")


    # OptimizableVariable around ScaledPolynomial
    const OVS = OptimizableVariable{N, DT, P,  B, V} where {N, DT, P<: ScaledPolynomial,  B, V}
    
    
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
                count_violations(PolynomialWrappers.coeffs(_ub),
                 PolynomialWrappers.coeffs(_lb), lb_violation_fun) > 0 && error("Lower boundary $(_lb) is higher than $(_ub)") 
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

    coeffs(o::OVS) = PolynomialWrappers.coeffs(o.p)

    lb_coeffs(ov::OVS)  =  PolynomialWrappers.coeffs(ov.lb)
    ub_coeffs(ov::OVS)  =  PolynomialWrappers.coeffs(ov.ub)

    derivative!(ov_der::OVS, ov::OVS) = PolynomialWrappers.derivative!(ov_der.p, ov.p)

    # const POSSIBLE_TAGS = (:lam, :C, )
"""
    Type to to store the single inverse problem (simplest case of the problem )
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
                                        initial_distribution::Union{OptimizableVariable, Number},
                                        thermocouples_locations::AbstractVector{DT},
                                        C,
                                        λ, 
                                        dλdT,
                                        thickness::Number, 
                                        xpoints_number::Int, 
                                        time_points_number::Union{Int,Nothing} = nothing;
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
            temperatures_needed = 3 - is_upper_flux_provided - is_lower_flux_provided
            issorted(time_data) || error("Time data must be sorted in ascending order")
            issorted(thermocouples_locations) || error("Thermocouple locations must be sorted in ascending order")
            NT = length(thermocouples_locations) # number of couples points
            all(Base.Fix2(<=, thickness), thermocouples_locations) || error("Thermocouple locations should be smaller than the value of thickness")
            NT < temperatures_needed && error("There should be at least $(temperatures_needed) thermocouples to solve the inverse problem")
            length(time_data) == size(temperatures, 1) || error("Number of rows in temperature data should be the same as the numbe rof time points")
            
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

            TN = temperatures_needed
            N = t_points
            ProblemType = typeof(direct_problem)
            DV = typeof(Tdata_evaluated)

            # all possibly optimizable variables are arranged into named tuple 

            #optimizable =(;λ = λ , C = C, dλdT = dλdT) # q_up = upper_flux , q_dwn = lower_flux , T₀ = initial_distribution, )
            
            optimizable =(; (k => v for (k, v) in zip((:λ,:C,:q_up, :q_dwn,:T₀, :dλdT), (λ, C, upper_flux, lower_flux, initial_distribution, dλdT)) if isa(v, OptimizableVariable))...)
            index_mapper = create_index_mapping(optimizable)
            O = typeof(optimizable)
            ON = length(optimizable)
            IM = typeof(index_mapper)
            new{DT, TN, N , ProblemType , CV , RG , DV, O, ON , IM}(        
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

        end
    end
    """
    fill_residual!(p::SingleInverseProblem)

Function  solves the direct problem and refills the resiaduals matrix
"""
function fill_residual!(p::SingleInverseProblem)
        @. p.residual = p.Tdata_evaluated - p.Tdata_measured
        return nothing
end

solve_direct_problem!(p::SingleInverseProblem) = solve_problem!(p.direct_problem)
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
    """
    update_all_optimizables!(p::SingleInverseProblem, x_vector::AbstractVector)

Function updates all optimizable variables with respect to the flags vectors
"""
function update_all_optimizables!(p::SingleInverseProblem, x_vector::AbstractVector)
        _modify_optimizables(p.optimizable, p.index_mapper , x_vector)
        is_λ_optimizable(p) && modify_λ_derivative!(p)
        return nothing
    end
@unroll function _modify_optimizables(optimizables , index_mapper , x)
 @unroll for i  in 1 : length(index_mapper)
            (ov, r) = (optimizables[i] , index_mapper[i])
            #modify!(ov, view(x , r))
            modify!(ov, x , r)
         end
end
    optimizable_parnumber(p::SingleInverseProblem) = sum(optimizable_parnumber, p.optimizable)

    is_λ_optimizable(p::SingleInverseProblem) = haskey(p.optimizable,:λ)

    function modify_λ_derivative!(p::SingleInverseProblem) 
        derivative!(p.optimizable.dλdT, p.optimizable.λ)
    end

    """
    fill_starting_vector(p::SingleInverseProblem{DT}) where DT

Scans all optimizable variables and returns the tuple with `(;x₀ - starting vector, lb - lower boundaries vector, ub - upper boundaries vector)``  for the optimization
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
function regularization_loss( p::SingleInverseProblem{DT,TN,N,P,CV,RG,  DV ,  O, ON}) where {DT , TN , N , P,CV, RG <: FiniteDifferenceRegularization,  DV ,  O, ON}
        return sum(finite_difference_regularization_loss , p.optimizable)
    end

include("covariances.jl")

    """
    discrepancy(x , p::SingleInverseProblem{DT}) where DT

Evaluates the weighted least-sqaure discrepancy of the corresponding inverse problem 
"""
function discrepancy(x , p::SingleInverseProblem{DT}) where DT

        update_all_optimizables!(p , x) # refreshes the values of parameters without solving the direct problem 
        solve_direct_problem!(p) # solves the direct problem 
        fill_residual!(p) # fills residual matrix 
        loss = covariance_loss(p) # applies weighted least squares
        p.include_constraints_violation_to_loss[] && (loss += p.ψ[] * constraints_loss(p)) # adds constraints loss to the main discrepancy (if they are needed)
        loss += p.α[] * regularization_loss(p)
        return loss

    end

    function interpolate_matrix!(Mout, t , M , tnew , start_col::Int = 1, stop_col::Int = 0)
        
        stop_col <= 0 && (stop_col = size(M,2))
        iter_step = start_col <= stop_col ? 1 : -1
        for (i , c) in enumerate(eachcol(M)[start_col : iter_step : stop_col])
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
            if abs(locations_vector[counter + 1] - zero_shift - xi) <= atol
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
    const ALL_REGULARIZATION_TYPES = subtypes(AbstractRegularization)
    const ALL_COVARIANCE_TYPES = subtypes(AbstractCovariance)
end
