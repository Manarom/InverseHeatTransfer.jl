

import .OneDHeatTransfer as OHT
import .PolynomialWrappers as PW

function evaluate_jacobian_lam(p::HeatTransferProblem{D, CF,LF,LDF,ITF, G}) where {D,
                                                                                CF <: ScaledPolynomial,
                                                                                LF <: ScaledPolynomial,
                                                                                LDF,
                                                                                ITF,
                                                                                G <: AbstractGrid{N,M}} where {N,M}
    
    
    
    interpolators = OHT.compute_derivatives(p)
    dCdT = PW.derivative(p.C_f.fun)
    dLdT = PW.derivative(p.L_f.fun)
    
end