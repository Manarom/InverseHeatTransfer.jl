# covariances loss functions 
    covariance_loss(p::SingleInverseProblem{DT, TN, N,
                        PT , CV  } ) where {DT, TN, N,
                                            PT , CV <: NoCovariance } = sum(abs2 , p.residual)/(N * TN)
    """
        Covariance with diagonal elements provided externally as a function dependent on temperature 
    
    Callable object σ² returns the square of dispertion as a function of measured temperature, 
    cache - stores values of sigma once evaluated at the starting 
    """
    struct TemperatureDependentDiagonalCovariance{ DT , F, TN , N} <: AbstractCovariance 
        σ²::F
        cache::Matrix{DT}
        function TemperatureDependentDiagonalCovariance(σ² , Tij::AbstractMatrix{DT}) where DT
            cache = similar(Tij)
            (N , TN) = size(cache)
            @inbounds for j in 1:TN
                    @simd for i in 1:N
                        cache[i , j] = σ²(Tij[i , j])
                    end
            end
            return new{ DT , typeof(σ²), TN , N}(σ² , cache)
        end
    end 
    struct FixedDiagonalCovariance{DT , M , N} <: AbstractCovariance
        σ²::M
        FixedDiagonalCovariance(σ²::M) where {M <: AbstractVector{DT}} where {DT} = new{DT , M , length(σ²)}(σ²)
    end
    function covariance_loss(p::SingleInverseProblem{DT, TN, N,
                        PT , CV  } ) where {DT, TN, N,
                                            PT , CV <: FixedDiagonalCovariance{DT} }
        loss = zero(DT)
    
        # iteration over thermocouples
        @inbounds for j in 1 : TN
        # iteration over time
            @simd for i in 1 : N
                r = p.residual[i , j]
                loss += (r * r) / p.covariance.σ²[i]
            end
        end
        return loss/(N * TN)
    end


    struct AR1Covariance{DT} <: AbstractCovariance
        τ::DT
        σ²::DT
    end
    """
            covariance_loss(p::SingleInverseProblem{DT, TN, N,
                            PT , CV  } ) where {DT , 
                                                PT , CV <: AR1Covariance{DT}} where { TN, N}

            Evaluates covariance for the first order autoregression , [Ornstein - Uhlenbeck process](https://en.wikipedia.org/wiki/Ornstein%E2%80%93Uhlenbeck_process)

        Measured temperature is assumed coorrelated with the following covariance:
        Σᵢⱼ = exp( - (i - j)*Δt/τ) 
        where `Δt` is the time step of the uniform grid (thus is works only on uniform grid solvers)
        `τ` - is the relaxation time 
        In this case there is an analytical inversion formula:
        `rᵗΣ⁻¹r = [r₁² + rₙ² + (1 + ρ²)Σ|i=2:n-1|(rᵢ²) - 2ρΣ|i=1:n-1|(rᵢrᵢ₊₁)]/(1 - ρ²)`
        `ρ = exp(- Δt/τ)`
    """
    function covariance_loss(p::SingleInverseProblem{DT, TN, N,
                        PT , CV  } ) where {DT , 
                                            PT <: HeatTransferProblem , 
                                            CV <: AR1Covariance{DT}} where {TN, N}

        !isa(p.direct_problem.grid , UniformGrid) && error("Grid must be uniform")
        dt = timestep(p.direct_problem)
        # as far as this type of covariance takes into account the timestep value need to recalculate 
        ρ = exp(- dt / p.covariance.τ)
        inv_den = 1.0 / (p.covariance.σ² * (1.0 - ρ^2))
        ρ2_plus_1 = 1.0 + ρ^2
        two_ρ = 2.0 * ρ   
        loss = zero(DT)

        @inbounds for j in StaticInt(1) : StaticInt(TN)

            s_squares = zero(DT) # squared amplitudes
            s_cross = zero(DT) # cross products
            
            @inbounds for i in 2 : N - 1
                r_curr = p.residual[i, j]
                s_squares += r_curr * r_curr # rᵢ²
                s_cross += r_curr * p.residual[i + 1 , j] # rᵢrᵢ₊₁
            end

            s_squares *= ρ2_plus_1 #  (1 + ρ²)Σ|i=2:n-1|(rᵢ²)
            # need to add cross product for i = 1 than mult by two_ρ
            @inbounds begin
                r1 = p.residual[1, j]
                rN = p.residual[N, j]
                r2 = p.residual[2, j]
                s_cross += r1 * r2
                s_squares += r1^2 + rN^2
            end
            s_cross *=  two_ρ
            loss += (s_squares  - s_cross) * inv_den
           # @show loss
        end
    
        return loss/(N * TN)
   end
    """
    covariance_loss(::SingleInverseProblem)

Function applies associated coavriance to and evaluate  the loss 
see `ALL_COVARIANCE_TYPES` for the list of implemented coavriances
"""
function covariance_loss(::SingleInverseProblem{DT, TN, N,
                                PT , CV  } ) where {DT, TN, N,
                                PT , CV  }
                                 error("Covariance $(CV) is not implemented") 
    end
    """
    fill_covariance_cache!(::SingleInverseProblem)

By default covariance fill cache do nothing
"""
    function fill_covariance_cache!(::SingleInverseProblem) end
    function fill_covariance_cache!(::SingleInverseProblem{DT, TN, N, PT , CV } ) where {DT, 
                                                        PT , 
                                                        CV <: FixedDiagonalCovariance{ DT , F , NCOV} } where {F, TN, N , NCOV} 
        @assert NCOV == N "Covariance weighting vector length $(NCOV) and residual column vector length $(N)" 
    end
    function fill_covariance_cache!(p::SingleInverseProblem{DT, TN, N,
                                    PT , CV  } ) where {DT, 
                                                        PT , 
                                                        CV <: TemperatureDependentDiagonalCovariance{ DT , F , TNCOV , NCOV} } where {F, TN, N , TNCOV , NCOV}
            @assert NCOV == N "Covariance weighting vector length $(NCOV) and residual column vector length $(N)" 
            @assert TNCOV == TN "Covariance weighting vector length $(TNCOV) and residual column vector length $(TN)" 
            @inbounds for j in 1:TN
            # iteration over time
                @simd for i in 1:N
                    Tij = p.Tdata_measured[i , j]
                    p.covariance.cache[i , j] = p.covariance.σ²(Tij)
                end
            end      

    end
    function covariance_loss(p::SingleInverseProblem{DT, TN, N,
                        PT , CV  } ) where {DT, 
                                            PT , CV <: TemperatureDependentDiagonalCovariance{F , DT , TN , N} } where {F, TN, N}
        
        loss = zero(DT)
        # iteration over thermocouples
        @inbounds for j in 1 : TN
            @simd for i in 1 : N
                r = p.residual[i , j]
                sigma_sq = p.covariance.cache[i , j]
                loss += (r * r) / sigma_sq
            end
        end
        return loss/(N * TN)
    end
    """
        Covariance proportional to the value of temperature to take into account relative 
    accuracy of temperature measurements
    """
    struct RelativeDiagonalCovariance{DT} <: AbstractCovariance
        relative_sigma::DT # value
        floor_sigma::DT    # 
    end
    function covariance_loss(p::SingleInverseProblem{DT, TN, N,
                        PT , CV  } ) where {DT, TN, N,
                                            PT , CV <: RelativeDiagonalCovariance{DT} }
        loss = zero(DT)
        
        rel_s = p.covariance.relative_sigma
        floor_s = p.covariance.floor_sigma
    
        # iteration over thermocouples
        @inbounds for j in 1 : TN
        # iteration over time
            @simd for i in 1 : N
                Tij = p.Tdata_measured[i , j]
                sigma_sq = (rel_s * Tij)^2 + floor_s^2
                r = p.residual[i , j]
                loss += (r * r) / sigma_sq
            end
        end
        return loss/(N * TN)
    end


    # trying to implement the lazy iterator over  residuals
    """ResidualIterator{W , PT , CV, N , TN , DT} structure to iterate over residauls as a vector without allocating 
    new vector 
    W - Val{true} - weighted residuals, unweighted otherwise
    PT - problem type 
    CV - porblem covariance type 
    N - number of time steps 
    TN - number of residual vector columns 
    DT - data type 
    RT - residual matrix eltype 
    """
    struct ResidualIterator{W , PT , CV, N , TN , DT , RT}
        p::PT
        r::RT
       function  ResidualIterator(  ::W, 
                                    p::PT) where PT <: SingleInverseProblem{DT, TN, N,
                                    DP , CV  } where { W <: Union{Val{true} , Val{false}},
                                                       DT, TN, N,
                                                       DP , CV  }

            new{W , PT , CV , N , TN , DT , typeof(p.residual)}(p , p.residual)
            
        end
    end

    Base.length(::ResidualIterator{W , PT , CV, N , TN}) where {W , PT , CV, N , TN} = N * TN
    Base.eltype(::ResidualIterator{W , PT , CV, N , TN , DT}) where {W , PT , CV, N , TN , DT} = DT

    function Base.iterate(iter::ResidualIterator{W , PT , CV, N }, state=1) where {W , PT , CV, N }
        if state > length(iter)
            return nothing
        end 
        i = (state - 1) % N + 1
        j = (state - 1) ÷ N + 1
        return (_get_weighted_residual_val(iter , i , j) , state + 1)
    end
    _get_weighted_residual_val(iter::ResidualIterator{Val{false}}, i , j) = iter.r[i, j]
    _get_weighted_residual_val(iter::ResidualIterator{Val{true}, PT , CV} , i , j) where {PT , CV <: NoCovariance} = iter.r[i, j]
    _get_weighted_residual_val(iter::ResidualIterator{Val{true}, PT , CV} , i , j) where {PT , CV <:TemperatureDependentDiagonalCovariance}  = iter.r[i , j] / sqrt(iter.p.covariance.cache[i , j])
    function _get_weighted_residual_val(iter::ResidualIterator{Val{true}, PT , CV} , i , j) where {PT , CV <: AR1Covariance{DT}} where DT 
        dt = timestep(iter.p.direct_problem)
        ρ = exp(- dt / iter.p.covariance.τ)
        return if i == 1
            iter.r[1 , j] / sqrt(iter.p.covariance.σ²)
        else
            σ_eps = sqrt(iter.p.covariance.σ² * (one(DT) - ρ^2))
            (iter.r[i , j] - ρ * iter.r[i - 1 , j]) / σ_eps
        end
    end
    function _get_weighted_residual_val(iter::ResidualIterator{Val{true}, PT , CV} , i , j) where {PT , CV <: RelativeDiagonalCovariance{DT}} where DT
        rel_s = iter.p.covariance.relative_sigma
        floor_s = iter.p.covariance.floor_sigma
        sigma_sq = (rel_s * Tij)^2 + floor_s^2
        return iter.r[i , j] / sqrt(sigma_sq)
    end

    struct ResidualColumn{W, PT, CV, N, TN, DT, RT} <: AbstractVector{DT} # {W , PT , CV, N , TN , DT , RT}
        iter::ResidualIterator{W, PT, CV, N, TN, DT, RT} # iterator over weighted residuals 
        col_idx::Int 
    end

    Base.size(::ResidualColumn{W, PT, CV, N}) where {W, PT, CV, N} = (N,)
    Base.IndexStyle(::Type{<:ResidualColumn}) = IndexLinear()

    @inline function Base.getindex(rc::ResidualColumn{W, PT, CV, N}, i::Int) where {W, PT, CV, N}
        global_state = i + (rc.col_idx - 1) * N
        return first(Base.iterate(rc.iter, global_state))
    end

    struct ResidualCols{W, PT, CV, N, TN, DT, RT}
        iter::ResidualIterator{W, PT, CV, N, TN, DT, RT}
    end

    Base.length(::ResidualCols{W, PT, CV, N, TN}) where {W, PT, CV, N, TN} = TN
    Base.eltype(::ResidualCols{W, PT, CV, N, TN, DT, RT}) where {W, PT, CV, N, TN, DT, RT} = ResidualColumn{W, PT, CV, N, TN, DT, RT}

    function Base.iterate(rcols::ResidualCols, state_col=1)
        if state_col > length(rcols)
            return nothing
        end
        return (ResidualColumn(rcols.iter, state_col), state_col + 1)
    end

    Base.eachcol(iter::ResidualIterator) = ResidualCols(iter)