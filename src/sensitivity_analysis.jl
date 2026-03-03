

using .OneDHeatTransfer
using .PolynomialWrappers

function evaluate_jacobian_lam(p::HeatTransferProblem{D, CF,LF,LDF,ITF, G}) where {D,
                                                                                CF <: ScaledPolynomial,
                                                                                LF <: ScaledPolynomial,
                                                                                LDF,
                                                                                ITF,
                                                                                G <: AbstractGrid{N,M}} where {N,M}
    (Tx_mat, Txx_mat) = compute_derivatives(p)
    T_mat = 
    itp_T   = LinearInterpolation(nodes, T_mat)
    itp_Tx  = LinearInterpolation(nodes, Tx_mat)
    itp_Txx = LinearInterpolation(nodes, Txx_mat)
    
    dCdT = PolynomialWrappers.derivative(p.C_f.fun)
    dLdT = PolynomialWrappers.derivative(p.L_f.fun)
    
end