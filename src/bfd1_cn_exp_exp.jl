
#=
            fill_LHS!(LHS, Fm1, F, Fp1, m, solver_scheme, problem)
            fill_RHS!(b, D, Tm, Tmm1, Fm1, F, Fp1, phi, solver_scheme, problem)
            apply_bc!(LHS, b, bc_fun,  F, Tm, phi ,  m,  solver_scheme, problem)
            ldiv!(Tmp1,LHS,b) # solving
=#
"""
    fill_LHS!(LHS, Fm1, F, Fp1, _ , ::BFD1_CN_EXP_EXP, _ )

LHS for `m + 1`'th time step
`[... , -0.5*F, (1 + F),  - 0.5*F , ...]T⁺ `
"""
function fill_LHS!(LHS, Fm1, F, Fp1, _ , ::BFD1_CN_EXP_EXP, _ ) 
    fill_tridiag!(LHS, Fm1, F, Fp1, 1.0, -0.5, 1.0, -0.5)
end
"""
    fill_RHS!(b, D, Tm, _ , _ , F, _ , phi, ::BFD1_CN_EXP_EXP, _ )


`MT +  F*ϕ*(Tn+1 - Tn-1)² `
M - tridiagonal with stencil
M = `[... , 0.5*F, (1 - F),  0.5*F , ...]T `
"""
function fill_RHS!(b, D, Tm, _ , _ , F, _ , phi, ::BFD1_CN_EXP_EXP, _ ) 
            mul!(b,D,Tm) 
            @. b = b^2
            @. b *= F*phi
            # b =b + M*c as M = Tridiagonal(F*ap1, a0 + a*F, F*ap1)
            # column_sym_tridiag_muladd!(b::AbstractVector{T},c::AbstractVector{T},
            #                        F::AbstractVector{T}, a0, a, a1)
            # @. b += Tm # Tm + \vec{b}
            a0 = 1.0
            a = -1.0
            a1 = 0.5
            column_sym_tridiag_muladd!(b, Tm, F, a0, a, a1) # b = b + M*Tm , M = [..., a1*F, a0 + a*F , a1*F ,...]
            return nothing
end


#-----------------------------------DirichletBC---UPPER-------------------------------------------
"""
    apply_bc!(LHS, b,  bc_fun::BoundaryFunction{D, BF, V, <: DirichletBC , <: UpperBC}, _ , _ , _ , m::Int , ::BFD1_CN_EXP_EXP, problem::HeatTransferProblem)   where {D , BF , V}

`apply_bc!(LHS, b, bc_fun,  F, Tm, phi ,  m,  solver_scheme, problem)`

``[1, 0, ...] T⁺ = T⁺``

"""
function apply_bc!(LHS, b,  bc_fun::BoundaryFunction{D, BF, V, <: DirichletBC , <: UpperBC}, _ , _ , _ , m::Int , ::BFD1_CN_EXP_EXP, problem::HeatTransferProblem)   where {D , BF , V}
        LHS.d[1] = 1.0
        LHS.du[1] = 0.0
        b[1] = bc_fun(tvalue(problem,m + 1)) # 1st order BC upper, evaluating bc for Tm+1
        return nothing
end
#------------------------------------DirichletBC--LOWER-------------------------------------------
"""
    apply_bc!(LHS, b,  bc_fun::BoundaryFunction{D, BF, V, <: DirichletBC , <: LowerBC}, _ , _ ,_,  m::Int, ::BFD1_CN_EXP_EXP, problem::HeatTransferProblem)   where {D , BF , V}

    `apply_bc!(LHS, b, bc_fun,  F, Tm, phi ,  m,  solver_scheme, problem)`
`` [..., 0 , 1] T⁺ = T⁺ ``
"""
function apply_bc!(LHS, b,  bc_fun::BoundaryFunction{D, BF, V, <: DirichletBC , <: LowerBC}, _ , _ ,_,  m::Int, ::BFD1_CN_EXP_EXP, problem::HeatTransferProblem)   where {D , BF , V}
        LHS.d[end] = 1.0
        LHS.dl[end] = 0.0
        #b[end] = bc_fun(m + 1) # 1st order BC upper, evaluating bc for Tm+1
        b[end] = bc_fun(tvalue(problem,m + 1)) # 1st order BC upper, evaluating bc for Tm+1
        return nothing
end



#--------------------------------NeumanBC---UPPER----------------------------------------------
"""
    apply_bc!(LHS, b,  bc_fun::BoundaryFunction{D, BF, V, <: NeumanBC , <: UpperBC}, F , Tm ,  m::Int, ::BFD1_CN_EXP_EXP, problem::HeatTransferProblem)   where {D,F,V}

`apply_bc!(LHS, b, bc_fun,  F, Tm, phi ,  m,  solver_scheme, problem)`

``[(1 +  F) , -F , ...] T⁺ =  [(1 - F), F,...]*T - F*Δx*(q⁺ + q)/λ + 4F*ϕ*(Δx*q/λ)² ``

"""
function apply_bc!(LHS, b,  bc_fun::BoundaryFunction{D, BF, V, <: NeumanBC , <: UpperBC}, F , Tm , phi,  m::Int, ::BFD1_CN_EXP_EXP, problem::HeatTransferProblem)   where {D , BF , V}
            f = F[1]
            ϕ = phi[1]
            Ti = Tm[1]
            dx_div_λ = xstep(problem) / thermal_conductivity(problem, Ti)
            

            LHS.du[1] = - f
            LHS.d[1] = 1.0 + f

            
            #qp1 = bc_fun(m + 1)
            #q = bc_fun(m)
            tmp1 = tvalue(problem, m + 1)
            tm = tvalue(problem, m )
            qp1 = bc_fun(tmp1)
            q = bc_fun(tm)

            b[1] = (1 - f) * Tm[1] + f * Tm[2] -  f * dx_div_λ * (qp1 + q)+ 4.0 * f * ϕ * (dx_div_λ * q)^2.0 
            return nothing
end
#----------------------------------NeumanBC----LOWER-------------------------------------------
"""
    apply_bc!(LHS, b,  bc_fun::BoundaryFunction{D, BF, V, <: NeumanBC , <: LowerBC}, F , Tm , phi,  m::Int, ::BFD1_CN_EXP_EXP, problem::HeatTransferProblem)   where {D,F,V}

`apply_bc!(LHS, b, bc_fun,  F, Tm, phi ,  m,  solver_scheme, problem)`

``[... , -F , (1 + F)] T⁺ =  [..., F, (1 - F)]*T - F*Δx*(q⁺ + q)/λ + 4F*ϕ*(Δx*q/λ)²``

"""
function apply_bc!(LHS, b,  bc_fun::BoundaryFunction{D, BF, V, <: NeumanBC , <: LowerBC}, F , Tm , phi,  m::Int, ::BFD1_CN_EXP_EXP, problem::HeatTransferProblem)   where {D , BF , V}
    
            f = F[end]
            ϕ = phi[end]
            Ti = Tm[end]
            dx_div_λ = xstep(problem) / thermal_conductivity(problem, Ti)
            

            LHS.dl[end] = -  f
            LHS.d[end] = 1.0 + f

            
            #qp1 = bc_fun(m + 1)
            #q = bc_fun(m)
            tmp1 = tvalue(problem, m + 1) # next time
            tm = tvalue(problem, m ) # current time
            qp1 = bc_fun(tmp1)
            q = bc_fun(tm)
            b[end] = f*Tm[end-1] + (1 - f) * Tm[end] - f * dx_div_λ * (qp1 + q) + 4.0 * f * ϕ * (dx_div_λ * q)^2.0 

            return nothing
end

#--------------------------------RobinBC---UPPER----------------------------------------------
"""
    apply_bc!(LHS, b,  bc_fun::BoundaryFunction{D, BF, V, <: RobinBC , <: UpperBC}, F , Tm ,  m::Int, ::BFD1_CN_EXP_EXP, problem::HeatTransferProblem)   where {D,F,V}

    `apply_bc!(LHS, b, bc_fun,  F, Tm, phi ,  m,  solver_scheme, problem)`

    ```[(1 +  2F) , -2F , 0, ...] T⁺ =  T - 2F*Δx*q⁺/λ + 4F*ϕ*(Δx*q/λ)² ```

"""
function apply_bc!(LHS, b,  bc_fun::BoundaryFunction{D, BF, V, <: RobinBC , <: UpperBC}, F , Tm , phi, 
                 m::Int, ::BFD1_CN_EXP_EXP, problem::HeatTransferProblem)   where {D , BF , V}
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
    apply_bc!(LHS, b,  bc_fun::BoundaryFunction{D, BF, V, <: NeumanBC , <: LowerBC}, F , Tm , phi,  m::Int, ::BFD1_CN_EXP_EXP, problem::HeatTransferProblem)   where {D,F,V}

    `apply_bc!(LHS, b, bc_fun,  F, Tm, phi ,  m,  solver_scheme, problem)`

    ```[..., 0 , -2F , (1 +  2F)] T⁺ =  T - 2F*Δx*q⁺/λ + 4F*ϕ*(Δx*q/λ)² ```

"""
function apply_bc!(LHS, b,  bc_fun::BoundaryFunction{D, BF, V, <: RobinBC , <: LowerBC}, F , Tm , phi,  m::Int, ::BFD1_CN_EXP_EXP, problem::HeatTransferProblem)   where {D , BF , V}
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