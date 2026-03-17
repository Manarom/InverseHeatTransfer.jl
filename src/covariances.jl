# covariances loss functions 
    covariance_loss(p::SingleInverseProblem{DT, TN, N,
                        PT , CV  } ) where {DT, TN, N,
                                            PT , CV <: NoCovariance } = sum(abs2 , p.residual)/(N * TN)
    """
        Covariance with diagonal elements provided externally as a function dependent on temperature 
    
    Callable object σ² returns the square of dispertion as a function of measured temperature, 
    cache - stores values of sigma once evaluated at the starting 
    """
    struct TemperatureDependentDiagonalCovariance{F , DT , TN , N} <: AbstractCovariance 
        σ²::F
        cache::Matrix{DT}
        function TemperatureDependentDiagonalCovariance(σ² , Tij::AbstractMatrix{DT}) where DT
            cache = similar(Tij)
            (N , TN) = size(cache)
            @inbounds for j in 1:TN
                    @simd for i in 1:N
                        cache[i , j] = σ²(Tij[i,j])
                    end
            end
            return new{typeof(σ²) , DT , TN , N}(σ² , cache)
        end
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
ρ = exp(- Δt/τ)


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

        #@show inv_den

        @inbounds for j in StaticInt(1) : StaticInt(TN)
           # @show j 

            s_squares = zero(DT) # squared amplitudes
            s_cross = zero(DT) # cross products
            
            @inbounds for i in 2:N-1
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

By default covariance fill cachedo nothing
"""
function fill_covariance_cache!(::SingleInverseProblem) end
function fill_covariance_cache!(p::SingleInverseProblem{DT, TN, N,
                                PT , CV  } ) where {DT, 
                                                    PT , 
                                                    CV <: TemperatureDependentDiagonalCovariance{F , DT , TN , N} } where {F, TN, N}
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
        axes(p.Tdata_measured) == axes(cov.cache) || throw(DimensionMismatch("Dimentions of covariance cache and T_measured are somehow different"))
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
