# runner for julia codes and benchmarking
using BenchmarkTools,Plots,Polynomials
 #plotly()
gr()
using AllocCheck#, Revise
include("OneDHeatTransfer.jl")

lam_pars = [0.44, 0.21e-2, 0.35e-4]
lam_fun = Polynomials.ImmutablePolynomial(lam_pars) #T -> lam_poly(T) # теплопроводность
lam_der_poly = Polynomials.ImmutablePolynomial( length(lam_pars)> 1 ? derivative(lam_fun) : [0.0])
lam_der = Polynomials.ImmutablePolynomial( lam_der_poly)#T->lam_der_poly(T)#;% производная теплопроводности
plot(range(200.0,1000,30),lam_fun.(range(200.0,1000.0,30)))
plot(range(200.0,1000,30),lam_der.(range(200.0,1000.0,30)))
#plot(linspace(200,1000,30),lam_der(linspace(200,1000,30)));title("Производная теплопроводности, Вт/(м*К^2)")
SHOW_BENCHMARKS = true
N = 50#;% число точек сетки по координате
M = 50000#;% число точек сетки по времени
M_explicit = 50000
Tmax = 1000.0#; % максимальная температура
tmax = 100.0#; % режим нагрева
Tinit = 20.0#; % начальная температура
Cp = 1000.0#; % теплоемкость
Ro = 2700.0#;% плотность
H = 15e-3#; % толщина слоя в мм
@eval Cp_fun(_)= $Cp*$Ro#;% не зависит от температуры
@eval initT_f(_) =  $Tinit #;% стартовая температура постоянна
BC_dwn_f = Polynomials.ImmutablePolynomial([Tinit])#;% температура снизу постоянна
BC_up_f =  Polynomials.ImmutablePolynomial([Tinit, (Tmax - Tinit)/tmax, 1e-2])  #;% температура сверху линейно возрастает

#plot(linspace(0,tmax,100),BC_up_f(linspace(0,tmax,100)))##;title("Режим нагрева")
#% решаем диффур

u_BC_type = OneDHeatTransfer.DirichletBC()
l_BC_type = OneDHeatTransfer.DirichletBC()


(T,) = OneDHeatTransfer.BFD1_exp_exp_exp(Cp_fun, lam_fun,lam_der, 
                        H, tmax,initT_f,
                        BC_up_f,BC_dwn_f,
                        M_explicit,
                        N, 
                        upper_bc_type = u_BC_type, lower_bc_type = l_BC_type)


plot(T,st=:surface)


p_exp = OneDHeatTransfer.HeatTransferProblem(Cp_fun,lam_fun,lam_der,initT_f, 
                            H, N, tmax, M,
                            BC_up_f, BC_dwn_f, Float64,
                            u_BC_type,
                            l_BC_type )

s_exp = OneDHeatTransfer.BFD1_EXP_EXP_EXP()
    
OneDHeatTransfer.unified_fd_solver!(p_exp,s_exp)
plot(p_exp.T, st=:surface)
plot(p_exp.T .- T,st=:surface)


#(T2,g,bc_up,bc_dwn) = OneDHeatTransfer.BFD1_imp_exp_exp(Cp_fun, lam_fun,lam_der, H, tmax,initT_f,BC_up_f,BC_dwn_f,M,N)

#implicit solver 

p_imp = OneDHeatTransfer.HeatTransferProblem(Cp_fun,lam_fun,lam_der,initT_f, 
                            H, N, tmax, M,
                            BC_up_f, BC_dwn_f, Float64,
                            u_BC_type,
                            l_BC_type )

s_imp = OneDHeatTransfer.BFD1_IMP_EXP_EXP()
    
OneDHeatTransfer.unified_fd_solver!(p_imp,s_imp)
plot(p_imp.T, st=:surface)
plot(p_imp.T .- T,st=:surface)
# testing crank nicolson
p_cn = OneDHeatTransfer.HeatTransferProblem(Cp_fun,lam_fun,lam_der,initT_f, 
                            H, N, tmax, M,
                            BC_up_f, BC_dwn_f, Float64,
                            u_BC_type,
                            l_BC_type )

s_cn = OneDHeatTransfer.BFD1_CN_EXP_EXP()

#FDSolverScheme(OneDHeatTransfer.BFD1(),
 #                                  OneDHeatTransfer.CN(),
 #                                  OneDHeatTransfer.EXP_NL(),
 #                                  OneDHeatTransfer.EXP_NL())
    
OneDHeatTransfer.unified_fd_solver!(p_cn,s_cn)
plot(p_cn.T,st=:surface)
plot(p_cn.T .- p_imp.T,st=:surface)

# testing BFD2_IMP_EXP_EXP
p_bfd2_imp = OneDHeatTransfer.HeatTransferProblem(Cp_fun,lam_fun,lam_der,initT_f, 
                            H, N, tmax, M,
                            BC_up_f, BC_dwn_f, Float64,
                            u_BC_type,
                            l_BC_type )

s_bfd2_imp = OneDHeatTransfer.BFD2_IMP_EXP_EXP()                            
OneDHeatTransfer.unified_fd_solver!(p_bfd2_imp,s_bfd2_imp)
plot(p_bfd2_imp.T,st=:surface)
plot(p_bfd2_imp.T .- p_cn.T,st=:surface)

# testing BFD2_CN_EXP_EXP
p_bfd2_cn = OneDHeatTransfer.HeatTransferProblem(Cp_fun,lam_fun,lam_der,initT_f, 
                            H, N, tmax, M,
                            BC_up_f, BC_dwn_f, Float64,
                            u_BC_type,
                            l_BC_type )
                            
s_bfd2_cn = OneDHeatTransfer.BFD2_CN_EXP_EXP()                            
OneDHeatTransfer.unified_fd_solver!(p_bfd2_cn,s_bfd2_cn)
plot(p_bfd2_cn.T, st=:surface)
plot(p_bfd2_cn.T .- p_bfd2_imp.T,st=:surface)

if SHOW_BENCHMARKS
    bm1 = @benchmark OneDHeatTransfer.BFD1_exp_exp_exp($Cp_fun, $lam_fun,$lam_der, H, tmax,$initT_f,$BC_up_f,$BC_dwn_f,M,N)
    bm_exp = @benchmark OneDHeatTransfer.unified_fd_solver!($p_exp,$s_exp)
    bm_imp = @benchmark OneDHeatTransfer.unified_fd_solver!($p_imp,$s_imp)
    bm_cn = @benchmark OneDHeatTransfer.unified_fd_solver!($p_cn,$s_cn)
    bm_bfd2_imp = @benchmark OneDHeatTransfer.unified_fd_solver!($p_bfd2_imp,$s_bfd2_imp)
    bm_bfd2_cn = @benchmark OneDHeatTransfer.unified_fd_solver!($p_bfd2_cn,$s_bfd2_cn)

    println("Direct solver (explicit scheme) NxM= $(N)x$(M)")
    display(bm1)
    println("Explicit scheme NxM= $(N)x$(M)")
    display(bm_exp)
    println("Implicit NxM= $(N)x$(M)")
    display(bm_imp)
    println("Crank-nicolson NxM= $(N)x$(M)")
    display(bm_cn)
    println("BFD2-implicit NxM= $(N)x$(M)")
    display(bm_bfd2_imp)
    println("BFD2 - Crank-Nicolson NxM= $(N)x$(M)")
    display(bm_bfd2_cn)
end
#BFD1_imp_exp_exp(Cp_fun, lam_fun,lam_der, H, tmax,initT_f,BC_up_f,BC_dwn_f,M,N)

