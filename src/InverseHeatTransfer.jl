module InverseHeatTransfer
    using LinearAlgebra,Reexport
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
    struct SingleInverseProblem{DT <: Number, TN , N ,
                    ProblemType <: HeatTransferProblem ,
                    CV <: AbstractCovariance, 
                    RG<: AbstractRegularization } # TN - couples number, N - timesteps number
        
        thermocouple_locations::SVector{TN , DT} # coordinates of all thermocouples
        thermocouple_indices::SVector{TN - 2 , DT} # indices of internal thermocouples in problem TMAT  - temperature distribution matrix 
        thermocouple_values # values of measured temperatures over time , number of rows  - N, 
        # number of columns must be equal to the number of locations 
        total_thickness
        direct_problem::ProblemType
        covariance::CV # covariance matrix
        regularization::RG # regularization matrix
        residual_vector # raw residual vector
        function SingleInverseProblem()
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
    function loss_function(u, ip::PropertyInversion{TP, N, T})
        fill_params!( u, ip)
        for (i, p_i) in enumerate(ip.problems)
            

        end
    end
end
