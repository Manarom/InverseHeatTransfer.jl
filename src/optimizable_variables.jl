
    abstract type AbstractRegularization end
    struct NoRegularization <: AbstractRegularization end
    struct FiniteDifferenceRegularization <: AbstractRegularization end
    struct FixedDiagonalRegularization <: AbstractRegularization end
    """
            Wrappes any modifiable and callcable variable of type P, which has parameters accessabel by `coeffs` function
    and can be bounded by `lb` and `ub` constraints, constraints violation can be checked  with  `lb_violation_fun`
    and `ub_violation_fun` functions. The parameters of `ub` and `lb` objects should be accessable by `lb_coeffs` and 
    `ub_coeffs`, there also should be a function to evaluate the derivative of `lb` and `ub` with respect to their 
    parameters.
    
        Currently implemented for `ScaledPolynomial` from `ScaledPolynomials`
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

        Wrapper interface for some callable object of type `P `which can be mutated by index `flag`. This
        structure can be bounded from the `top` and from the `bottom` using some objects `ub` and `lb` , 
        both of the same type `V`. Variables `iu` and `il` are just flags which can be used to turn `on` 
        and `off` this boundaries, boundaries are compared to the current state of the varible `p` using
        `f1` and `f2` functions, which must take two argument of `P` and `V` types and return a single value.

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
    """
        diagonal_regularization_loss(ov::OV{N,DT})

    Simple regularization 
    """
    function fixed_diagonal_regularization_loss(ov::OV{N,DT} ) where {N , DT}

        loss = zero(DT)
        avg_coeffs = zero(DT)
        @inbounds @simd for i in 1 : N
            _c = coeffs(ov)[i]
            #avg_coeffs += _c
            loss +=  _c^2.0
        end 
        
        return loss / N 
    end
    fview_coeffs(ov::OV) = view(coeffs(ov), ov.flag)
    fview_lb_coeffs(ov::OV) = view(lb_coeffs(ov), ov.flag)
    fview_ub_coeffs(ov::OV) = view(ub_coeffs(ov), ov.flag)
    fview_lb_violation(ov::OV)  = view(ov.lb_violation, ov.flag)
    fview_ub_violation(ov::OV)  = view(ov.ub_violation, ov.flag)

    extract_all_params(ov::OV) = copy(coeffs(ov))
    extract_optimizable_params(ov::OV) = copy(fview_coeffs(ov))
    # interface
    # necessary function for details see the implementation for ScaledPolynomials
    coeffs(::OV{N,DT,P}) where {N,DT,P}  = error("OptimizableVariable wrappers is not implemenented for $(P) type")
    lb_coeffs(::OV{N,DT,P})  where {N,DT,P}  = error("OptimizableVariable wrappers around $(P) is not implemenented")
    ub_coeffs(::OV{N,DT,P}) where {N,DT,P}  = error("OptimizableVariable wrappers around $(P) is not implemenented")
    # not necessary
    derivative!(::OV, ::OV)  = error("OptimizableVariable wrappers around is not implemented")