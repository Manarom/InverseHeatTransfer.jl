### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 9bb428c9-5947-4bc2-b4bc-d05959850af2
begin
import Pkg
	Pkg.activate(@__DIR__)
	Pkg.instantiate()
	using PlutoUI,Interpolations, BenchmarkTools,PlutoPlotly,Interpolations,StaticArrays,RecipesBase,Plots, Revise, PrettyTables
	import  InverseHeatTransfer.ScaledPolynomials as PW
end

# ╔═╡ c7400800-f152-11f0-a386-29e067f9216c


# ╔═╡ 5fe54406-20d7-4747-bcd4-c5c0f17af47a


# ╔═╡ 9386589f-1767-46a8-8a60-ec3f44a0970d
solver_path = joinpath(@__DIR__,"..","src","solvers")

# ╔═╡ a91ae0d5-6a90-472e-abd7-a5b02cf470fe
polynomials_path = joinpath(@__DIR__,"..","src","polynomials")

# ╔═╡ ada6d5ce-ac5e-4de8-b9fe-4394be6a34ca
includet(joinpath(solver_path,"OneDHeatTransfer.jl"))

# ╔═╡ 1936a00f-69f3-48e4-9176-6b2186e40c2c
	begin 
		#PW = Main.ScaledPolynomials
		OHT = Main.OneDHeatTransfer
	end

# ╔═╡ 6c87fcb7-55e6-44f4-8978-12e74b9cfdbf
begin
	#N = 50# coordinate grid
	#M = 5000# time grid
	Tmax = 2000.0# max temperature
	tmax = 100.0# heating time
	Tinit = 20.0# starting temperature
	Cp = 1000.0# heat capacity
	Ro = 2700.0# density
	H = 15e-3# layer thickness in m
	PolyType = PW.BernsteinSymPoly
	#@eval Cp_fun(_) = $Cp*$Ro# capacity
	#@eval initT_f(x) = $Tinit*cos(0*x*pi/$H) # starting 
end;

# ╔═╡ d6145328-5f43-4148-8911-232c9fed39c4
md"""
	Select solver: $(@bind solver_name  Select(collect(Main.OneDHeatTransfer.AVAILABLE_SCHEMES)))

	"""

# ╔═╡ d8e97e88-da66-48d7-b5cd-4c1b3d9a4c4c
begin 
	fd_solver =OHT.unified_fd_solver!
	solver_type_type = OHT.AVAILABLE_SCHEMES[solver_name]
	solver_type = solver_type_type()
end

# ╔═╡ 57b23767-880e-4375-b9b4-bd0164df4a4b
md"""

	Select upper BC type: $(@bind upper_bc_type  Select(subtypes(Main.OneDHeatTransfer.AbstractBoundaryCondition)))
	"""

# ╔═╡ eeb3d1a0-be52-444d-aa82-89c639d5aece
begin 
	is_upper_dirichle =  upper_bc_type == OHT.DirichletBC
	is_upper_neuman =  upper_bc_type == OHT.NeumanBC  
	ubc_mult = !is_upper_dirichle ? 5e6 : Tmax
	def_val_ubc = !is_upper_dirichle  ? 1e5/2 : 1e3
end;

# ╔═╡ 36e2c7b2-0cb0-4684-8d01-03e107349cb2
@bind  upper_bc_pars PlutoUI.combine() do Child
	md"""
	upper BC(T) = a₀x³ + a₁x² + a₂x + a₃
	
	a0 = $(
		Child(Slider(ubc_mult*(-1:0.01:1), default = 20.0, show_value = true))
	) \
	a1 = $(
		Child(Slider(ubc_mult*(-1:0.01:1), default = def_val_ubc, show_value = true))
	)\
	a2 = $(
		Child(Slider(ubc_mult*(-1:0.01:1), default = def_val_ubc, show_value = true))
	)\
	a3 = $(
		Child(Slider(ubc_mult*(-1:0.01:1), default = def_val_ubc, show_value = true))
	)\
	a0 = $(
		Child(Slider(ubc_mult*(-1:0.01:1), default = 20.0, show_value = true))
	) \
	a1 = $(
		Child(Slider(ubc_mult*(-1:0.01:1), default = def_val_ubc, show_value = true))
	)\
	a2 = $(
		Child(Slider(ubc_mult*(-1:0.01:1), default = def_val_ubc, show_value = true))
	)\
	a3 = $(
		Child(Slider(ubc_mult*(-1:0.01:1), default = def_val_ubc, show_value = true))
	)\
	"""
end

# ╔═╡ 9db7fc50-8d3e-416b-8a57-2c819d6dc426
md""" 
Plot ditribution over $(@bind plot_type Select([:time,:coordinate], default = :coordinate))

Abscissa max value $(@bind x_max_value Slider(0.1:0.1:tmax, default = tmax, show_value = true) )

"""

# ╔═╡ f9f5058d-6943-4d96-9821-eafd89100fc6
md"""

	Select lower BC type: $(@bind lower_bc_type  Select(subtypes(Main.OneDHeatTransfer.AbstractBoundaryCondition),default = subtypes(Main.OneDHeatTransfer.AbstractBoundaryCondition)[2]))

	"""

# ╔═╡ 4d348c80-1035-4561-85e2-02c9bf12486a
begin 
	is_lower_dirichle =  lower_bc_type ==OHT.DirichletBC
	is_lower_neuman =  lower_bc_type == OHT.NeumanBC
	lbc_mult = !is_lower_dirichle ? 5e5 : 1e3
	def_val_lbc = !is_lower_dirichle ? 0 : 0
end;

# ╔═╡ 8c0b27e0-7f0b-4921-a679-ceda2ad46e82
@bind  lower_bc_pars PlutoUI.combine() do Child
	md"""
	lower  BC(T) = a₀x³ + a₁x² + a₂x + a₃
	
	a0 = $(
		Child(Slider(lbc_mult*(-1:0.01:1), default = def_val_lbc, show_value = true))
	) \
	a1 = $(
		Child(Slider(lbc_mult*(-1:0.01:1), default = def_val_lbc, show_value = true))
	)\
	a2 = $(
		Child(Slider(lbc_mult*(-1:0.01:1), default = def_val_lbc, show_value = true))
	)\
	a3 = $(
		Child(Slider(lbc_mult*(-1:0.01:1), default = def_val_lbc, show_value = true))
	)\
	"""
end

# ╔═╡ e746f979-e39d-4838-8ae1-a64a1afed1f4
@bind  initT_pars PlutoUI.combine() do Child
	md"""
	initial temperature distribution
	``T(x) = a_0x^3 + a_1x^2 + a_2x + a_3``
	
	a0 = $(
		Child(Slider(-100:1.0:2000,default=20.0,show_value = true))
	) \
	a1 = $(
		Child(Slider(-100:1.0:2000,default=20.0,show_value = true))
	)\
	a2 = $(
		Child(Slider(-100:1.0:2000,default=20.0,show_value = true))
	)\
	a3 = $(
		Child(Slider(-100:1.0:2000,default=20.0,show_value = true))
	)\
	"""
end

# ╔═╡ ddc9ba47-c4f2-4c77-8a9c-c598b7eb322b
@bind  c_pars PlutoUI.combine() do Child
	md"""
	``Cp(T) = a_0T^3 + a_1T^2 + a_2T + a_3``
	
	a0 = $(
		Child(Slider(100:1:2000,default=600,show_value = true))
	) \
	a1 = $(
		Child(Slider(100:1:2000,default=1100,show_value = true))
	)\
	a2 = $(
		Child(Slider(100:1:2000,default=1200,show_value = true))
	)\
	a3 = $(
		Child(Slider(100:1:2000,default=1200,show_value = true))
	)\
	"""
end

# ╔═╡ 4316370d-583d-4b48-8632-d17bc34209bd
md"""
	time nods = $(@bind M Slider(2:2:1000, default = 50, show_value = true))

	coordinate nod number =  $(@bind N Slider(2:2:500, default = 50, show_value = true))
	"""

# ╔═╡ ea882e26-5972-443f-a7ac-3075701b90fe
begin 
	N_t_dummy_points = 100
	t_dummy_range = range(0,tmax, N_t_dummy_points)
	t_dummy = collect(t_dummy_range);
	T_dummy = collect(range(Tinit,Tmax,N))
	x_dummy = collect(range(0,H,N))
end;

# ╔═╡ da943c51-3f33-4f10-bc3d-9129f76399df
begin 

	is_upper_time_dependent = is_upper_dirichle || is_upper_neuman
	(ubc_small_name,up_x_data_name, up_x_data) = 
		if is_upper_time_dependent 
			("t","Time,s", t_dummy) 
		else 
			("T","Temperature, ᵒC", T_dummy)
		end
	up_x_data = is_upper_time_dependent ? t_dummy : T_dummy
	(min1,max1) =  extrema(up_x_data)
	BC_up_f = PW.ScaledPolynomial(PolyType, ntuple( _-> 0.0, length(upper_bc_pars)), min1 , max1)

	
end;

# ╔═╡ 1611afbd-f13c-4532-aec9-6404b33e3f15
begin 

	is_lower_time_dependent = is_lower_dirichle || is_lower_neuman
	(lbc_small_name,low_x_data_name, low_x_data) = 
		if is_lower_time_dependent 
			("t","Time,s", t_dummy) 
		else 
			("T","Temperature, ᵒC", T_dummy)
		end
	low_x_data = is_lower_time_dependent ? t_dummy : T_dummy
	BC_dwn_f = PW.ScaledPolynomial(PolyType,(0.0,0.0,0.0,0.0), extrema(low_x_data)...)

end;

# ╔═╡ a4cd1c4f-f744-4ec1-bd4d-52a9e83a73ff
begin 
	initT_f = PW.ScaledPolynomial(PolyType, initT_pars, extrema(x_dummy)...)
	init_values = initT_f.(x_dummy)
	Plots.plot(1e3*x_dummy,init_values,title="Starting temperature",label= nothing,linewidth = 4)
	xlabel!("Coordinate, mm")
	ylabel!("Temperature distribution, ᵒC")
end

# ╔═╡ 0eae0fd4-8b5d-448a-b935-e06d469f080b
@bind  lam_pars PlutoUI.combine() do Child
	md"""
	``\lambda(T) = a_0T^3 + a_1T^2 + a_2T + a_3``
	
	a0 = $(
		Child(Slider(1e-3:0.01:100,default=17,show_value = true))
	) \
	a1 = $(
		Child(Slider(1e-3:0.01:100,default=12,show_value = true))
	)\
	a2 = $(
		Child(Slider(1e-3:0.01:100,default=8,show_value = true))
	)\
	a3 = $(
		Child(Slider(1e-3:0.01:100,default=10,show_value = true))
	)\
	"""
end

# ╔═╡ 84bfb85d-def3-4a08-be60-01de2d68be36
md" Use plotly backend for surf $(@bind is_use_plotly CheckBox(default = false))"

# ╔═╡ 8992ffda-fcbb-4994-9acb-1080056903b7
md"""
	camera angle α $(@bind α Slider(-180:180,default = 0))
	camera angle β $(@bind β Slider(-180:180,default = 0))
	"""

# ╔═╡ 331ee2df-931b-4267-8795-1fe402c7fdaa
@bind T_range PlutoUI.combine() do Child 
md"""
	``T_{min} \ ^oC``= $(
		Child(Slider(-200:1.0:1000,default=-200,show_value = true))
	) \
	``T_{max} \ ^oC`` = $(
		Child(Slider(0:1.0:3000,default=2000,show_value = true))
	)\
	"""
end

# ╔═╡ b7cab05e-c24f-4268-ba7b-f361c1e9076d
begin 
	lam_T_range = range(T_range...,length = 100)
	lam_fun = PW.ScaledPolynomial(PolyType,lam_pars, T_range...)
end;

# ╔═╡ ed45f038-ae2b-4888-8547-ab9638170842
begin 
	Cp_fun = PW.ScaledPolynomial(PolyType,ntuple( i->Ro*c_pars[i], 4), T_range...)
	cp_values = Cp_fun.(lam_T_range)
	p_cp = Plots.plot(lam_T_range,cp_values./Ro,title="Specific heat",label= nothing,linewidth = 4)
	xlabel!(p_cp,"Temperature, ᵒC")
	ylabel!(p_cp,"Specific heat , J/(kg*K)")
end

# ╔═╡ 8f6d3be5-ffb5-41d5-8e50-f4353193c53c
begin 
	lam_values = lam_fun.(lam_T_range)
	p_lam = Plots.plot(lam_T_range,lam_values,title="thermal conductivity",label= nothing,linewidth = 4)
	xlabel!(p_lam,"Temperature, ᵒC")
	ylabel!(p_lam,"Thermal conductivity, W/(m*K)")
end

# ╔═╡ 38976547-9399-49f1-b83e-0833327ea0ed
begin 
	lam_der = PW.derivative(lam_fun)
problem = Main.OneDHeatTransfer.HeatTransferProblem(Cp_fun, 					lam_fun,lam_der,initT_f, H,N, tmax, M, BC_up_f, BC_dwn_f, Float64, upper_bc_type(),lower_bc_type());

end

# ╔═╡ 9fbe7c0d-90b7-4348-8b22-d8058452a02b
begin 
	
	lower_bc_pars
	upper_bc_pars
	PW.refill!(BC_up_f, upper_bc_pars)
	PW.refill!(BC_dwn_f,lower_bc_pars)
	
    fd_solver(problem,solver_type)

	# lower BC plot
	p_lbc = Plots.plot(low_x_data,BC_dwn_f.(low_x_data),title="lower BC:"*string(lower_bc_type),label= nothing,linewidth = 4)
	xlabel!(p_lbc,low_x_data_name)
	ylabel!(p_lbc,"lower BC")
	# upper BC plot
	p_ubc = Plots.plot(up_x_data,BC_up_f.(up_x_data),title="upper BC:"*string(upper_bc_type),label= nothing,linewidth = 4)
	xlabel!(p_ubc,up_x_data_name)
	ylabel!(p_ubc,"upper BC")

	
	TCN = problem.T 
end;

# ╔═╡ 3a88786c-cefd-4daa-bc71-1d7cf04b9a2e
p_ubc

# ╔═╡ a554174d-82d6-4acd-89a9-bc21aca78e7d
p_lbc

# ╔═╡ ec6e4678-a807-47f6-845d-dc529f38e80c
begin 
	NN = 9
	step  = M ÷ NN
	coord = Main.OneDHeatTransfer.xrange(problem)
	time =  Main.OneDHeatTransfer.trange(problem)
	
	Tinterp = linear_interpolation((coord, time), TCN)
	bknd = is_use_plotly ?  Plots : Plots
	is_plot_coordinate = plot_type == :coordinate
	x_iterpolation = is_plot_coordinate ? range(0.0,minimum((tmax,x_max_value) ), length = 11) : range(0.0,minimum((H,x_max_value*1e-3)), length = 11)

	if is_plot_coordinate
		ppp = bknd.plot(grid = true,gridlinewidth=3,gridstyle = :dot,minorgrid=true,legend_position = :best,legend_title = "time:")
		for t_cur in x_iterpolation
			bknd.plot!(ppp,1e3*coord,Tinterp.(coord,t_cur), label = string(round(t_cur,digits = 2)),linewidth = 3,linealpha = 0.8)
		end
		xlabel!(ppp,"Coordinate, mm")
		ylabel!(ppp,"Temperature")
	else # plot distribution over time 
		ppp = bknd.plot(grid = true,gridlinewidth=3,gridstyle = :dot,minorgrid=true,legend_position = :best,legend_title = "coord:")
		for x_cur in x_iterpolation
			bknd.plot!(ppp,time,Tinterp.(x_cur,time), label = string(round(1e3*x_cur,digits = 2)),linewidth = 3,linealpha = 0.8)
		end		
		xlabel!(ppp,"Time, s")
		ylabel!(ppp,"Temperature")
	end
	surf_view = @view TCN[:,1:step:M]
	ppp
end

# ╔═╡ 3acc0be7-dabc-4577-917d-26530cf192bd
if is_use_plotly 
	#plotly()
	tr = PlutoPlotly.surface(z=TCN, colorscale="Heat")
	layout = Layout(
    width=800, 
    height=600, 
    autosize=true,
    margin=attr(l=0, r=0, b=0, t=50),  # Minimize margins
    scene=attr(
        camera=attr(eye=attr(x=1.5, y=1.5, z=1.2)),  
        aspectmode="cube"  
    )
)
	 PlutoPlotly.plot(tr)
else
	#gr()
	Plots.plot(surf_view,st = :surface,camera=(α, β),grid = true,framestyle = :box, cmap=:hot)
end

# ╔═╡ e4aac0a1-d9c0-4e20-8f86-151090a9651f
begin 
	lam_pars
	lam_prs = PW.BernsteinSymPoly(lam_pars)
	#PW.polyfit!(lam_prs,lam_T_range,lam_values)
	#lam_der = PW.derivative( lam_prs)#;% производная теплопроводности
end;

# ╔═╡ 7a22e820-13df-4e56-8507-d1b089fe7f52
md" Use bench? $(@bind bench CheckBox(default = false))"

# ╔═╡ 0ed065e8-4f97-46ba-bf07-eb02fb2ae9d6
!bench || @benchmark fd_solver($problem,$solver_type)

# ╔═╡ 5e7504fd-29f1-41a3-b940-5921c471a4f5
problem

# ╔═╡ 802f203a-8f5b-4a64-a26b-4fe7f812ff56
solver_type

# ╔═╡ 3f7a81bf-6cd4-43b4-b5db-8b891f4e4722
md" Compare to reference? $(@bind reference_eval CheckBox(default = false))"

# ╔═╡ 2742b8e1-4a9e-4021-87a2-c11ab15503cd
md"""
	camera angle α $(@bind αref Slider(-180:180,default = 45, show_value = true))

	camera angle β $(@bind βref Slider(-180:180,default = 45, show_value = true))
	"""

# ╔═╡ e2f4dde1-e256-4980-81ee-a61d582f1b47
md" 

Benchmark all solvers (takes much time ) $(@bind is_bench_all CheckBox(default = false))

Compare accuracy for all solvers $(@bind is_check_all CheckBox(default = false))
"

# ╔═╡ eae7841c-1b5b-497d-a67d-c7502f2d26ed
if reference_eval 
	fnames = [:mean_ΔT, :Tmax, :std_ΔT, :max_ΔT, :tmax,:xmax, :Fo]
	md" Select quantity to plot $(@bind quant2plot Select(fnames, default = :max_ΔT))"
end

# ╔═╡ 79a82e5d-417a-4168-96a6-aed2f2190a8c
md" Select benchmark output $(@bind bm_field Select([:time,:memory,:allocs,:gctime] ))"

# ╔═╡ 18e22e41-1b9d-4e3b-a12d-92f943e8b17a
md"""
	xscale = $(@bind bm_xscale Select([:identity, :log10]))
	yscale = $(@bind bm_yscale Select([:identity, :log10]))
	"""

# ╔═╡ 2cd30b6e-409a-4f5b-a53d-8ba3b7510b0a
md" Select benchmark function $(@bind bm_fun Select([mean => :mean,maximum => :maximum] ))"

# ╔═╡ 7017c716-8682-4a74-9bfe-7dc09879db16
function evaluate_stat(problem,pref_T_interpolation,solver_name)
	coord = OHT.xrange(problem)
	time =  OHT.trange(problem)
	Tref = similar(problem.T)
	@inbounds for (i,x) in enumerate(coord)
		for (j,t) in enumerate(time)
			Tref[i,j]= pref_T_interpolation(x,t)
		end
	end
	Tdiff = Tref .- problem.T	
	meanT = mean(Tdiff)
	stddT = std(Tdiff)
	maxdT = maximum(abs,Tdiff)
	out = findmax(abs,Tdiff)
	tmax = time[out[2][2]]
	xmax = 1e3*coord[out[2][1]]
	Tmax = problem.T[out[2]]
	Fo = OHT.fourier_number(problem,Tmax)
	stat_str = md"""
	for $(solver_name)  finite-difference scheme with 

	Δt = $(OHT.timestep(problem))

	Δx = $(1e3*OHT.xstep(problem))

	average ΔT ± std is  $(meanT) ± $(stddT) 

	max ΔT is $(maxdT) at 

	t = $(tmax) s 

	x = $(xmax) mm

	"""
	return (;stat_str = stat_str, Tref = Tref, ΔT = Tdiff, mean_ΔT = meanT, std_ΔT = stddT, max_ΔT = maxdT, tmax = tmax, xmax = xmax, name = solver_name, Tmax = Tmax, Fo = Fo)
end

# ╔═╡ e525b2ec-f4a7-4052-af2b-7a9b4fa4f520
function interpolate_pref(p_ref)
	x2int = OHT.xrange(p_ref)
	y2int = OHT.trange(p_ref)
	T_ref_int = cubic_spline_interpolation((x2int,y2int), p_ref.T)
end

# ╔═╡ 60506dec-faf8-4ce1-9d2f-d5ea816a7825
if reference_eval
	reference_problem = Main.OneDHeatTransfer.HeatTransferProblem(Cp_fun, 					lam_fun,lam_der,initT_f, H,500, tmax,50000,BC_up_f,BC_dwn_f,Float64,upper_bc_type(),lower_bc_type());

#=
problem = Main.OneDHeatTransfer.HeatTransferProblem(Cp_fun, 					lam_fun,lam_der,initT_f, H,N, tmax, M, BC_up_f, BC_dwn_f, Float64, upper_bc_type(),lower_bc_type());
	=#

	
	reference_solver_type = OHT.BFD1_CN_EXP_EXP()
	fd_solver(reference_problem, reference_solver_type)
	T_ref_int = interpolate_pref(reference_problem)
	stat_cur = evaluate_stat(problem, T_ref_int, solver_name)
	Tdiff = stat_cur.ΔT
	#=
	if is_use_plotly 
		tr2 = PlutoPlotly.surface(x = 1e3*coord, y = time, z=Tdiff, colorscale="Viridis")
		layout2 = Layout(
    		width=800, 
    		height=600, 
    		autosize=true,
    		margin=attr(l=0, r=0, b=0, t=50),  # Minimize margins
    		scene=attr(
        		xaxis_title="coordinate",
            	yaxis_title="time", 
            	zaxis_title="T - Tref"

    		)
		)
		p_res = PlutoPlotly.plot(tr2)
	else
		p_res = Plots.plot(1e3*coord, time,Tdiff, st = :surf, xlabel = "coordinate" , ylabel = "time", zlabel = "T - Tref",camera=(αref, βref),)
	end
	=#
end

# ╔═╡ 518d31b7-35d3-46e9-8ce4-fa634c1a036e
!reference_eval|| md"""$(stat_cur.stat_str)"""

# ╔═╡ 7ef0ec03-3c98-44d6-b237-1bf4deecf834
if reference_eval 
	Tref_bench = interpolate_pref(reference_problem)
	#m_vect = [10,20]
	m_vect = [5,10,50, 100,200, 300,600, 800,  1000, 2000, 5000, 10000]
	m_number = length(m_vect)
	schemes2compare = filter(t->t[1] != :BFD1_EXP_EXP_EXP, OHT.AVAILABLE_SCHEMES )
	sol_number = length(schemes2compare)
	row_number = length(m_vect)
	stats = Matrix{Any}(missing,(m_number, sol_number))
	sol_names = collect(keys(schemes2compare))
	accuracy_matrix = Matrix{Float64}(undef,(row_number, sol_number))
	bm_accuracy_matrix = Matrix{Float64}(undef,(row_number, sol_number))
	table_columns = vcat(:m,:Δt,sol_names)
end;

# ╔═╡ db4e5649-8724-4a6f-a12e-acc3fccc2f89
if is_check_all && reference_eval
	for i in 1:m_number
		mi = m_vect[i]
		p_i = OHT.HeatTransferProblem(Cp_fun,lam_fun,lam_der,initT_f, H,N, tmax,mi, BC_up_f,BC_dwn_f, Float64, upper_bc_type(), lower_bc_type())
		for (j,(name,type)) in enumerate(schemes2compare)
			OHT.unified_fd_solver!(p_i, type())
			bm_cur = !is_bench_all ? nothing : @benchmark OHT.unified_fd_solver!($p_i, $(type())) seconds = 2
			stats[i,j] = (evaluate_stat(p_i,Tref_bench, name), bm_cur)
		end
	end
end

# ╔═╡ d6c38053-4773-42e9-acd4-05ddb761c48a

if (is_check_all && reference_eval) 
	p_stat_M = Plots.plot(marker=true, markersize=8, markershape=:auto)
	quant2plot_vect = Vector{Float64}(undef,row_number)
	for j in 1 : sol_number # column  - solver type
		for i in 1:row_number# row - M numbers
			si = stats[i,j][1]
			quant2plot_vect[i] =  getfield(si,quant2plot)
		end
		accuracy_matrix[:,j] .= quant2plot_vect
		plot!(p_stat_M, m_vect, abs.(quant2plot_vect),yscale = :log10, xscale = :log10, label = string(sol_names[j]),marker=true, markersize=8, markershape= :auto, linewidth = 3)
	end
	xlabel!(p_stat_M, "Number of time steps")
	ylabel!(p_stat_M,string(quant2plot))
	title!(p_stat_M,"Maximum error in temperature estimation")
	p_stat_M	
end

# ╔═╡ c13976ff-0f0b-448c-a97e-575d9020985d
if (is_check_all && reference_eval)
	dt_vect = tmax./(m_vect .- 1)
		pretty_table(HTML,hcat(m_vect, dt_vect,accuracy_matrix),column_labels = table_columns, title = "$(quant2plot)")
end

# ╔═╡ 8a576b4f-5073-451f-b277-aa80284ad0ac
if is_check_all && reference_eval && is_bench_all
	bm_stat_M = Plots.plot(marker=true, markersize=8, markershape=:auto)
	bench2plot_vect = Vector{Float64}(undef,row_number)
	fun_on_bench = bm_fun
	field_of_bench = bm_field
	for j in 1 : size(stats,2) # column  - solver type
		for i in 1:row_number# row - M numbers
			bmi = stats[i,j][2] # benchmark obj
			bench2plot_vect[i] =  getfield(fun_on_bench(bmi),field_of_bench)
		end
		bm_accuracy_matrix[:,j] .= bench2plot_vect
		plot!(bm_stat_M, m_vect, bench2plot_vect, xscale = bm_xscale,yscale=bm_yscale, label = string(sol_names[j]),marker=true, markersize=8, markershape= :auto, linewidth = 3)
	end
	xlabel!(bm_stat_M, "Number of time steps")
	ylabel!(bm_stat_M,string(bm_fun)*"_"*string(bm_field)*" nanoseconds")
	title!(bm_stat_M,"Execution time")
	bm_stat_M	
end

# ╔═╡ 6f0a8a7b-4fcd-45ff-bccd-bae01ed53748
if is_check_all && reference_eval && is_bench_all
		pretty_table(HTML,hcat(m_vect, dt_vect,bm_accuracy_matrix*1e-6),column_labels = table_columns, title = "$(bm_fun) $(bm_field) ,ms")
end

# ╔═╡ 91a53ea2-1576-4651-a015-68b7ada03c9a
!(is_check_all && reference_eval && is_bench_all) || md"""
	Show bench for:

	solver : $(@bind show_bench_solver Select(sol_names))

	M= : $(@bind show_bench_m Select(m_vect))
	"""

# ╔═╡ 8833c701-d50b-4561-bda9-a93df065a179
if is_check_all && reference_eval && is_bench_all
	mind = findfirst(isequal(show_bench_m),m_vect)
	solind = findfirst(isequal(show_bench_solver), sol_names)
	stats[mind,solind][2]
end

# ╔═╡ d5762b64-6db0-4460-afa7-ac7e0b36e448
#=
pretty_table(HTML,hcat(["Direct","EmPoint", "AutoDiff"],discr_values,grad_values,hess_values),header = ["calculated using:","discrepancy", "gradient", "hessian"],top_left_str= "Table of discrepancy and its derivatives values, calculated for T=$(TtryK),K" )

=#

# ╔═╡ Cell order:
# ╠═9bb428c9-5947-4bc2-b4bc-d05959850af2
# ╠═c7400800-f152-11f0-a386-29e067f9216c
# ╠═5fe54406-20d7-4747-bcd4-c5c0f17af47a
# ╠═9386589f-1767-46a8-8a60-ec3f44a0970d
# ╠═a91ae0d5-6a90-472e-abd7-a5b02cf470fe
# ╠═ada6d5ce-ac5e-4de8-b9fe-4394be6a34ca
# ╟─1936a00f-69f3-48e4-9176-6b2186e40c2c
# ╠═6c87fcb7-55e6-44f4-8978-12e74b9cfdbf
# ╠═ea882e26-5972-443f-a7ac-3075701b90fe
# ╠═eeb3d1a0-be52-444d-aa82-89c639d5aece
# ╠═4d348c80-1035-4561-85e2-02c9bf12486a
# ╠═da943c51-3f33-4f10-bc3d-9129f76399df
# ╠═1611afbd-f13c-4532-aec9-6404b33e3f15
# ╠═d8e97e88-da66-48d7-b5cd-4c1b3d9a4c4c
# ╟─d6145328-5f43-4148-8911-232c9fed39c4
# ╟─57b23767-880e-4375-b9b4-bd0164df4a4b
# ╟─3a88786c-cefd-4daa-bc71-1d7cf04b9a2e
# ╟─36e2c7b2-0cb0-4684-8d01-03e107349cb2
# ╟─9db7fc50-8d3e-416b-8a57-2c819d6dc426
# ╟─f9f5058d-6943-4d96-9821-eafd89100fc6
# ╟─a554174d-82d6-4acd-89a9-bc21aca78e7d
# ╟─8c0b27e0-7f0b-4921-a679-ceda2ad46e82
# ╟─a4cd1c4f-f744-4ec1-bd4d-52a9e83a73ff
# ╟─e746f979-e39d-4838-8ae1-a64a1afed1f4
# ╟─ed45f038-ae2b-4888-8547-ab9638170842
# ╟─ddc9ba47-c4f2-4c77-8a9c-c598b7eb322b
# ╟─8f6d3be5-ffb5-41d5-8e50-f4353193c53c
# ╟─ec6e4678-a807-47f6-845d-dc529f38e80c
# ╟─4316370d-583d-4b48-8632-d17bc34209bd
# ╟─3acc0be7-dabc-4577-917d-26530cf192bd
# ╟─38976547-9399-49f1-b83e-0833327ea0ed
# ╟─0eae0fd4-8b5d-448a-b935-e06d469f080b
# ╟─9fbe7c0d-90b7-4348-8b22-d8058452a02b
# ╟─84bfb85d-def3-4a08-be60-01de2d68be36
# ╟─8992ffda-fcbb-4994-9acb-1080056903b7
# ╟─331ee2df-931b-4267-8795-1fe402c7fdaa
# ╟─b7cab05e-c24f-4268-ba7b-f361c1e9076d
# ╟─e4aac0a1-d9c0-4e20-8f86-151090a9651f
# ╟─7a22e820-13df-4e56-8507-d1b089fe7f52
# ╟─0ed065e8-4f97-46ba-bf07-eb02fb2ae9d6
# ╟─5e7504fd-29f1-41a3-b940-5921c471a4f5
# ╟─802f203a-8f5b-4a64-a26b-4fe7f812ff56
# ╟─3f7a81bf-6cd4-43b4-b5db-8b891f4e4722
# ╟─2742b8e1-4a9e-4021-87a2-c11ab15503cd
# ╟─60506dec-faf8-4ce1-9d2f-d5ea816a7825
# ╟─518d31b7-35d3-46e9-8ce4-fa634c1a036e
# ╟─e2f4dde1-e256-4980-81ee-a61d582f1b47
# ╟─7ef0ec03-3c98-44d6-b237-1bf4deecf834
# ╟─eae7841c-1b5b-497d-a67d-c7502f2d26ed
# ╟─db4e5649-8724-4a6f-a12e-acc3fccc2f89
# ╟─d6c38053-4773-42e9-acd4-05ddb761c48a
# ╟─c13976ff-0f0b-448c-a97e-575d9020985d
# ╟─79a82e5d-417a-4168-96a6-aed2f2190a8c
# ╟─18e22e41-1b9d-4e3b-a12d-92f943e8b17a
# ╟─2cd30b6e-409a-4f5b-a53d-8ba3b7510b0a
# ╟─8a576b4f-5073-451f-b277-aa80284ad0ac
# ╟─6f0a8a7b-4fcd-45ff-bccd-bae01ed53748
# ╟─91a53ea2-1576-4651-a015-68b7ada03c9a
# ╟─8833c701-d50b-4561-bda9-a93df065a179
# ╟─7017c716-8682-4a74-9bfe-7dc09879db16
# ╟─e525b2ec-f4a7-4052-af2b-7a9b4fa4f520
# ╟─d5762b64-6db0-4460-afa7-ac7e0b36e448
