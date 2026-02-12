#-----------------------------------DirichletBC---UPPER-------------------------------------------
"""
    apply_bc!(LHS, b,  bc_fun::BoundaryFunction{D, BF, V, <: DirichletBC , <: UpperBC}, _ , _ , _ , m::Int , ::AbstractFDScheme, problem::HeatTransferProblem)   where {D , BF , V}

    Dirichlet BC is are the same for all FD schemes

``[1, 0, ...] T⁺ = T⁺``

"""
function apply_bc!(LHS, b,  bc_fun::BoundaryFunction{D, BF, V, <: DirichletBC , <: UpperBC}, _ , _ , _ , m::Int , ::AbstractFDScheme, problem::HeatTransferProblem)   where {D , BF , V}
        LHS.d[1] = 1.0
        LHS.du[1] = 0.0
        b[1] = bc_fun(tvalue(problem,m + 1)) # 1st order BC upper, evaluating bc for Tm+1
        return nothing
end
#------------------------------------DirichletBC--LOWER-------------------------------------------
"""
    apply_bc!(LHS, b,  bc_fun::BoundaryFunction{D, BF, V, <: DirichletBC , <: LowerBC}, _ , _ ,_,  m::Int, ::AbstractFDScheme, problem::HeatTransferProblem)   where {D , BF , V}

    Dirichlet BC is are the same for all FD schemes
        
`` [..., 0 , 1] T⁺ = T⁺ ``
"""
function apply_bc!(LHS, b,  bc_fun::BoundaryFunction{D, BF, V, <: DirichletBC , <: LowerBC}, _ , _ ,_,  m::Int, ::AbstractFDScheme, problem::HeatTransferProblem)   where {D , BF , V}
        LHS.d[end] = 1.0
        LHS.dl[end] = 0.0
        #b[end] = bc_fun(m + 1) # 1st order BC upper, evaluating bc for Tm+1
        b[end] = bc_fun(tvalue(problem,m + 1)) # 1st order BC upper, evaluating bc for Tm+1
        return nothing
end
