
#=
            fill_LHS!(LHS, Fm1, F, Fp1, m, solver_scheme, problem)
            fill_RHS!(b, D, Tm, Tmm1, Fm1, F, Fp1, m, phi, solver_scheme, problem)
            apply_bc!(LHS, b, bc_fun,  F, Tm, phi ,  m,  solver_scheme, problem)
            ldiv!(Tmp1,LHS,b) # solving
=#
#const FIRST_STEP_BFD2 = BFD1_CN_EXP_EXP()
"""
LHS for `m + 1`'th time step
`[... , -0.5*F, (3/2 + F),  - 0.5*F , ...]T⁺ `
"""
function fill_LHS!(LHS, Fm1, F, Fp1, m , ::BFD2_CN_EXP_EXP, problem ) 
    m != 1 || return fill_LHS!(LHS, Fm1, F, Fp1, m , FIRST_STEP_BFD2, problem ) 
    fill_tridiag!(LHS, Fm1, F, Fp1, 1.5, -0.5 , 1.0 , -0.5)
    return nothing
end
"""
    fill_RHS!(b, D, Tm, Tmm1, Fm1, F, Fp1, m, phi, ::BFD2_CN_EXP_EXP, problem)
``0.5F*Tᵐn+1 +(2 - F)Tᵐn + 0.5F*Tᵐn-1 - 0.5*Tᵐ⁻¹ +  F*ϕ*(Tn+1 - Tn-1)² ``
"""
function fill_RHS!(b, D, Tm, Tmm1, Fm1 , F, Fp1 , m, phi, ::BFD2_CN_EXP_EXP, problem) 
            m != 1 || return fill_RHS!(b, D, Tm, Tmm1, Fm1, F, Fp1, m, phi, FIRST_STEP_BFD2, problem)
            mul!(b,D,Tm) 
            @. b = b^2
            @. b *= F*phi
            @. b -=  0.5 * Tmm1# Tm + \vec{b}
            # b = M*c + b
            # tridiag_muladd!(b, c , Fm1 ,
            #                        F , Fp1 , a0, 
            #                        am1, a, ap1)
            #tridiag_muladd!(b , Tm, Fm1, F, Fp1, 2.0, 0.5, -1.0, 0.5)
            #column_sym_tridiag_muladd!(b::AbstractVector{T},c::AbstractVector{T},
            #                        F::AbstractVector{T}, a0, a, a1)
            column_sym_tridiag_muladd!(b, Tm, F, 2.0, -1.0, 0.5)
            return nothing
end

#--------------------------------NeumanBC---UPPER----------------------------------------------
"""
    apply_bc!(LHS, b,  bc_fun::BoundaryFunction{D, BF, V, <: NeumanBC , <: UpperBC}, F , Tm ,  m::Int, ::BFD2_IMP_EXP_EXP, problem::HeatTransferProblem)   where {D,F,V}

    `apply_bc!(LHS, b, bc_fun,  F, Tm, phi ,  m,  solver_scheme, problem)`

    ```[(1.5 +  F) , -F ,  ...] T⁺ =  [(2 - F), F,... ]*T - 0.5Tm-1 - (FΔx/λ)*(q⁺ + q) + 4F*ϕ*(Δx*q/λ)² ```
"""
function apply_bc!(LHS, b,  bc_fun::BoundaryFunction{D, BF, V, <: NeumanBC , <: UpperBC}, F , Tm , phi,  m::Int, ::BFD2_CN_EXP_EXP, problem::HeatTransferProblem)   where {D , BF , V}
            
            m != 1 || return apply_bc!(LHS, b,  bc_fun, F , Tm , phi,  m, FIRST_STEP_BFD2 , problem) 

            f = F[1]
            ϕ = phi[1]
            Ti = Tm[1]
            dx_div_λ = xstep(problem) / thermal_conductivity(problem, Ti)
            

            LHS.du[1] = - f
            LHS.d[1] = 1.5 +  f

            
            #qp1 = bc_fun(m + 1)
            #q = bc_fun(m)
            tmp1 = tvalue(problem, m + 1)
            tm = tvalue(problem, m )
            qp1 = bc_fun(tmp1)
            q = bc_fun(tm)
            Tmm1 = problem.T[1, m - 1] # previous timestep temperature of the first node
            b[1] = (2.0  - f) * Tm[1] + f*Tm[2]- 0.5 * Tmm1 -  f * dx_div_λ * (qp1 + q) + 4.0 * f * ϕ * (dx_div_λ * q)^2.0 
            return nothing
end
#----------------------------------NeumanBC----LOWER-------------------------------------------
"""
    apply_bc!(LHS, b,  bc_fun::BoundaryFunction{D, BF, V, <: NeumanBC , <: LowerBC}, F , Tm , phi,  m::Int, ::BFD2_IMP_EXP_EXP, problem::HeatTransferProblem)   where {D,F,V}

    `apply_bc!(LHS, b, bc_fun,  F, Tm, phi ,  m,  solver_scheme, problem)`

    ```[..., 0 , - F , (1.5 +  F)] T⁺ =   [...,F, (2 - F)]*T - 0.5Tm-1 - FΔx/λ*(q⁺ + q) + 4F*ϕ*(Δx*q/λ)² ```

"""
function apply_bc!(LHS, b,  bc_fun::BoundaryFunction{D, BF, V, <: NeumanBC , <: LowerBC}, F , Tm , phi,  m::Int, ::BFD2_CN_EXP_EXP, problem::HeatTransferProblem)   where {D , BF , V}
            m != 1 || return apply_bc!(LHS, b,  bc_fun, F , Tm , phi,  m, FIRST_STEP_BFD2 , problem) 
            f = F[end]
            ϕ = phi[end]
            Ti = Tm[end]
            dx_div_λ = xstep(problem) / thermal_conductivity(problem, Ti)
            

            LHS.dl[end] = -  f
            LHS.d[end] = 1.5 +  f

            
            #qp1 = bc_fun(m + 1)
            #q = bc_fun(m)
            tmp1 = tvalue(problem, m + 1) # next time
            tm = tvalue(problem, m ) # current time
            qp1 = bc_fun(tmp1)
            q = bc_fun(tm)
            Tmm1 = problem.T[end, m - 1] # previous timestep temperature of the first node
            b[end] =f*Tm[end - 1] + (2.0 - f)*Tm[end] - 0.5*Tmm1 - f * dx_div_λ * (qp1 + q)  + 4.0 * f * ϕ * (dx_div_λ * q)^2.0 
            return nothing
end

#--------------------------------RobinBC---UPPER----------------------------------------------
"""
    apply_bc!(LHS, b,  bc_fun::BoundaryFunction{D, BF, V, <: RobinBC , <: UpperBC}, F , Tm ,  m::Int, ::BFD2_IMP_EXP_EXP, problem::HeatTransferProblem)   where {D,F,V}

    `apply_bc!(LHS, b, bc_fun,  F, Tm, phi ,  m,  solver_scheme, problem)`

    ```[(1 +  2F) , -2F , 0, ...] T⁺ =  T - 2F*Δx*q⁺/λ + 4F*ϕ*(Δx*q/λ)² ```

"""
function apply_bc!(LHS, b,  bc_fun::BoundaryFunction{D, BF, V, <: RobinBC , <: UpperBC}, F , Tm , phi, 
                 m::Int, ::BFD2_CN_EXP_EXP, problem::HeatTransferProblem)   where {D , BF , V}
            error("to do later")
            f = F[1]
            ϕ = phi[1]
            Ti = Tm[1]
            dx_div_λ = xstep(problem) / thermal_conductivity(problem, Ti)
            

            LHS.du[1] = - 2.0 * f
            LHS.d[1] = 1 + 2.0 * f

            
            qp1 = bc_fun(m + 1)
            q = bc_fun(m)
            b[1] = Tm[1] - 2.0 * f * dx_div_λ * qp1 + 4.0 * f * ϕ * (dx_div_λ * q)^2.0 
            return nothing
end
#--------------------------------RobinBC---LOWER----------------------------------------------
"""
    apply_bc!(LHS, b,  bc_fun::BoundaryFunction{D, BF, V, <: NeumanBC , <: LowerBC}, F , Tm , phi,  m::Int, ::BFD2_IMP_EXP_EXP, problem::HeatTransferProblem)   where {D,F,V}

    `apply_bc!(LHS, b, bc_fun,  F, Tm, phi ,  m,  solver_scheme, problem)`

    ```[..., 0 , -2F , (1 +  2F)] T⁺ =  T - 2F*Δx*q⁺/λ + 4F*ϕ*(Δx*q/λ)² ```

"""
function apply_bc!(LHS, b,  bc_fun::BoundaryFunction{D, BF, V, <: RobinBC , <: LowerBC}, F , Tm , phi,  m::Int, ::BFD2_CN_EXP_EXP, problem::HeatTransferProblem)   where {D , BF , V}
            error("to do later")

            f = F[end]
            ϕ = phi[end]
            Ti = Tm[end]
            dx_div_λ = xstep(problem) / thermal_conductivity(problem, Ti)
            

            LHS.dl[end] = - 2.0 * f
            LHS.d[end] = 1 + 2.0 * f

            
            qp1 = bc_fun(m + 1)
            q = bc_fun(m)
            b[end] = Tm[end] - 2.0 * f * dx_div_λ * qp1 + 4.0 * f * ϕ * (dx_div_λ * q)^2.0 

            return nothing
end