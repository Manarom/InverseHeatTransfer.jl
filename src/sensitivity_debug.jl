using Plots
using DifferentialEquations
using Polynomials
using Revise
plotly()
include(raw"./solvers/OneDHeatTransfer.jl")
include(raw"./polynomials/PolynomialWrappers.jl")
PW = PolynomialWrappers
OHT = OneDHeatTransfer
includet("sensitivity_analysis.jl")
T_range = (20.0,1000.0)
lam_pars = (0.44, 0.21, 0.35, 1.2)
lam_fun = PW.ScaledPolynomial(PW.BernsteinSymPoly, lam_pars, T_range...)
lam_der = PW.derivative(lam_fun)

Ro = 2700.0#;% плотность
Cp_fun= PW.ScaledPolynomial(PW.BernsteinSymPoly,(Ro*900,Ro*1100,Ro*1100,Ro*1100),T_range...)

plot(range(200.0,1000,30),lam_fun.(range(T_range...,30)))
plot(range(200.0,1000,30),lam_der.(range(T_range...,30)))
plot(range(200.0,1000,30),Cp_fun.(range(T_range...,30))./Ro)

N = 150#;% число точек сетки по координате
M = 2000#;% число точек сетки по времени
M_explicit = 5000
(Tinit,Tmax) = T_range
tmax = 100.0#; % режим нагрева

H = 15e-3#; % толщина слоя в мм

@eval initT_f(_) =  $Tinit #;% стартовая температура постоянна
BC_dwn_f = Polynomials.ImmutablePolynomial([Tinit])#;% температура снизу постоянна
BC_up_f =  Polynomials.ImmutablePolynomial([Tinit, (Tmax - Tinit)/tmax, 1e-2])  #;% температура сверху линейно возрастает

#plot(linspace(0,tmax,100),BC_up_f(linspace(0,tmax,100)))##;title("Режим нагрева")
#% решаем диффур

u_BC_type = OHT.DirichletBC()
l_BC_type = OHT.DirichletBC()



p = OHT.HeatTransferProblem(Cp_fun,lam_fun,lam_der,initT_f, 
                            H, N, tmax, M,
                            BC_up_f, BC_dwn_f, Float64,
                            u_BC_type,
                            l_BC_type )
   
OHT.solve_problem!(p)
T = OHT.temperature_field(p)
Tx = OHT.temperature_gradient(p)
x = OHT.xrange(p.grid)
t = OHT.trange(p.grid)
surface(x,t,T, st=:surface)
plot(Tx, st=:surface)
plot(p.T, st=:surface)

plot(Base.Fix1(T,0.0).(1:100))

