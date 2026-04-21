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

# ╔═╡ a4f32130-c129-426f-8092-a688e981aea8
begin 
	import Pkg
	Pkg.activate(@__DIR__)
	Pkg.instantiate()
	using Interpolations , StaticArrays , Plots  , Revise , PrettyTables
	using Optimization , OptimizationOptimJL , OptimizationNLopt , OptimizationMetaheuristics
	using LinearAlgebra , DelimitedFiles
	using PlutoUI , PlutoPlotly , OrderedCollections
	using HypertextLiteral
	using Revise
	import InverseHeatTransfer
	using DataFrames
	using CSV
	using RecipesBase
	using Dates

end

# ╔═╡ c7400800-f152-11f0-a386-29e067f9216c
#using Polynomials,Interpolations,Interpolations,Polynomials,LegendrePolynomials,StaticArrays,RecipesBase,Plots, FFTW

# ╔═╡ b6133376-a592-48ee-8d63-a8413e20088d
#using BenchmarkTools, ProfileCanvas, PlutoUI, PlutoPlotly

# ╔═╡ 5fe54406-20d7-4747-bcd4-c5c0f17af47a
#using Revise,AllocCheck, PrettyTables, Reexport, LinearAlgebra

# ╔═╡ e6216c49-09a8-45c5-8c3a-6d1f7281cc3e
#using Optimization, OptimizationOptimJL, Observables, FunctionWrappers

# ╔═╡ 79ff4f54-ce5d-4f30-81ea-0c6cbdadceb0
#using FiniteDiff, OrderedCollections, Tables

# ╔═╡ 8a53eac7-3fa8-4d33-bf86-3cfcbdf120c5
#using PlutoTables, Accessors, AccessorsExtra

# ╔═╡ c806aa46-a30e-4431-b2e1-c7131ce86615
source_path = joinpath(@__DIR__,"..","src")

# ╔═╡ 32ac05db-67f0-4638-9842-d5ef37442efa
#includet(joinpath(source_path, "InverseHeatTransfer.jl"))

# ╔═╡ 1936a00f-69f3-48e4-9176-6b2186e40c2c
begin 
	IHT = InverseHeatTransfer
	PW = IHT.ScaledPolynomials
	OHT = IHT.OneDHeatTransfer
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
	H = 7e-3# layer thickness in m
	PolyType = PW.BernsteinSymPoly
	#@eval Cp_fun(_) = $Cp*$Ro# capacity
	#@eval initT_f(x) = $Tinit*cos(0*x*pi/$H) # starting 
end;

# ╔═╡ 4316370d-583d-4b48-8632-d17bc34209bd
md"""
	time nods = $(@bind M Slider(2:2:10000, default = 5000, show_value = true))

	coordinate nod number =  $(@bind N Slider(2:2:1000, default = 100, show_value = true))
	"""

# ╔═╡ ea882e26-5972-443f-a7ac-3075701b90fe
begin 
	N_t_dummy_points = 100
	t_dummy_range = range(0,tmax, N_t_dummy_points)
	t_dummy = collect(t_dummy_range);
	T_dummy = collect(range(Tinit,Tmax,N))
	x_dummy = collect(range(0,H,N))
end;

# ╔═╡ d6145328-5f43-4148-8911-232c9fed39c4
md"""
	Select solver: $(@bind solver_name  Select(collect(OHT.AVAILABLE_SCHEMES)))

	"""

# ╔═╡ d8e97e88-da66-48d7-b5cd-4c1b3d9a4c4c
begin 
	fd_solver =OHT.unified_fd_solver!
	solver_type_type = OHT.AVAILABLE_SCHEMES[solver_name]
	solver_type = solver_type_type()
	#solver_type.OHT.	
		# = Main.OneDHeatTransfer.FDSolverScheme(Main.OneDHeatTransfer.BFD1(),
                                  #Main.OneDHeatTransfer.IMP(),
                                   #Main.OneDHeatTransfer.EXP_NL(),
                                   #Main.OneDHeatTransfer.EXP_NL()) =#
end

# ╔═╡ 57b23767-880e-4375-b9b4-bd0164df4a4b
md"""

	Select upper BC type: $(@bind upper_bc_type  Select(subtypes(OHT.AbstractBoundaryCondition)))
	"""

# ╔═╡ eeb3d1a0-be52-444d-aa82-89c639d5aece
begin 
	is_upper_dirichle =  upper_bc_type == OHT.DirichletBC
	is_upper_neuman =  upper_bc_type == OHT.NeumanBC  
	ubc_mult = !is_upper_dirichle ? 5e6 : Tmax
	def_val_ubc = !is_upper_dirichle  ? 1e5/2 : 1e3
end;

# ╔═╡ 36e2c7b2-0cb0-4684-8d01-03e107349cb2
@bind  upper_bc_pars confirm(PlutoUI.combine() do Child
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
		Child(Slider(ubc_mult*(-1:0.01:1), default = def_val_ubc, show_value = true))
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
end)

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
	BC_up_f = PW.ScaledPolynomial(PolyType, ntuple( _-> 0.0, length(upper_bc_pars)), min1,max1)

	
end;

# ╔═╡ 9db7fc50-8d3e-416b-8a57-2c819d6dc426
md""" 
Plot ditribution over $(@bind plot_type Select([:time,:coordinate], default = :coordinate))

Abscissa max value $(@bind x_max_value Slider(0.1:0.1:tmax, default = tmax, show_value = true) )

"""

# ╔═╡ f9f5058d-6943-4d96-9821-eafd89100fc6
md"""

	Select lower BC type: $(@bind lower_bc_type  Select(subtypes(OHT.AbstractBoundaryCondition),default = subtypes(OHT.AbstractBoundaryCondition)[2]))

	"""

# ╔═╡ 4d348c80-1035-4561-85e2-02c9bf12486a
begin 
	is_lower_dirichle =  lower_bc_type ==OHT.DirichletBC
	is_lower_neuman =  lower_bc_type == OHT.NeumanBC
	lbc_mult = !is_lower_dirichle ? 5e5 : 1e3
	def_val_lbc = !is_lower_dirichle ? 0 : 0
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

# ╔═╡ a4cd1c4f-f744-4ec1-bd4d-52a9e83a73ff
begin 
	initT_f = PW.ScaledPolynomial(PolyType, initT_pars, extrema(x_dummy)...)
	init_values = initT_f.(x_dummy)
	Plots.plot(1e3*x_dummy,init_values,title="Starting temperature",label= nothing,linewidth = 4)
	xlabel!("Coordinate, mm")
	ylabel!("Temperature distribution, ᵒC")
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

# ╔═╡ 0eae0fd4-8b5d-448a-b935-e06d469f080b
@bind  lam_pars PlutoUI.combine() do Child
	md"""
	``\lambda(T) = a_0T^3 + a_1T^2 + a_2T + a_3``
	
	a0 = $(
		Child(Slider(1e-3:0.01:30,default=2,show_value = true))
	) \
	a1 = $(
		Child(Slider(1e-3:0.01:30,default=2,show_value = true))
	)\
	a2 = $(
		Child(Slider(1e-3:0.01:30,default=2,show_value = true))
	)\
	a3 = $(
		Child(Slider(1e-3:0.01:30,default=2, show_value = true))
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
		Child(Slider(-200:1.0:1000,default=10,show_value = true))
	) \
	``T_{max} \ ^oC`` = $(
		Child(Slider(0:1.0:3000,default=1200,show_value = true))
	)\
	"""
end

# ╔═╡ 86a00c2e-4896-417d-93d0-bc69b5823957
begin 
	lam_T_range = range(T_range...,length = 100)
	lam_fun = PW.ScaledPolynomial(PolyType,lam_pars, T_range...)
	lam_der = PW.derivative(lam_fun)
end;

# ╔═╡ ed45f038-ae2b-4888-8547-ab9638170842
begin 
	Cp_fun = PW.ScaledPolynomial(PolyType,ntuple( i->Ro*c_pars[i], 4), T_range...)
	cp_values = Cp_fun.(lam_T_range)
	p_cp = Plots.plot(lam_T_range,cp_values./Ro,title="Specific heat",label= nothing,linewidth = 4)
	xlabel!(p_cp,"Temperature, ᵒC")
	ylabel!(p_cp,"Specific heat , J/(kg*K)")
end

# ╔═╡ 38976547-9399-49f1-b83e-0833327ea0ed
	problem = IHT.HeatTransferProblem(Cp_fun, 					lam_fun,lam_der,initT_f, H,N, tmax, M, BC_up_f, BC_dwn_f, Float64, upper_bc_type(),lower_bc_type());

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

# ╔═╡ ec6e4678-a807-47f6-845d-dc529f38e80c
begin 
	NN = 9
	step  = M ÷ NN
	coord = IHT.xrange(problem)
	time =  IHT.trange(problem)
	
	Tinterp = linear_interpolation((coord, time), TCN)
	bknd = is_use_plotly ?  Plots : Plots
	is_plot_coordinate = plot_type == :coordinate
	x_iterpolation = is_plot_coordinate ? range(0.0,minimum((tmax,x_max_value) ), length = 11) : range(0.0,minimum((H,x_max_value*1e-3)), length = 11)

	if is_plot_coordinate
		ppp = bknd.plot(grid = true,gridlinewidth=3,gridstyle = :dot,minorgrid=true,legend_position = :best,legend_title = "time:", box = :on)
		for t_cur in x_iterpolation
			bknd.plot!(ppp,1e3*coord,Tinterp.(coord,t_cur), label = string(round(t_cur,digits = 2)),linewidth = 3,linealpha = 0.8)
		end
		xlabel!(ppp,"Coordinate, mm")
		ylabel!(ppp,"Temperature")
	else # plot distribution over time 
		ppp = bknd.plot(grid = true,gridlinewidth=3,gridstyle = :dot,minorgrid=true,legend_position = :best,legend_title = "coord:", box = :on)
		for x_cur in x_iterpolation
			bknd.plot!(ppp,time,Tinterp.(x_cur,time), label = string(round(1e3*x_cur,digits = 2)),linewidth = 3,linealpha = 0.8)
		end		
		xlabel!(ppp,"Time, s")
		ylabel!(ppp,"Temperature")
	end
	surf_view = @view TCN[:,1:step:M]
	ppp
end

# ╔═╡ a554174d-82d6-4acd-89a9-bc21aca78e7d
p_lbc

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

# ╔═╡ 8f6d3be5-ffb5-41d5-8e50-f4353193c53c
begin 
	lam_values = lam_fun.(lam_T_range)
	p_lam = Plots.plot(lam_T_range,lam_values,title="thermal conductivity",label= nothing,linewidth = 4)
	xlabel!(p_lam,"Temperature, ᵒC")
	ylabel!(p_lam,"Thermal conductivity, W/(m*K)")
end

# ╔═╡ 95508f3b-73b5-4cc5-b6bd-bb4d949bf00c
md"""

## Now we are going to test the inverse problem and try to fit the obtained results

### Generating ``measured`` data 


"""

# ╔═╡ 1182c47f-3cfb-4aed-a597-c79a0a786d3c
possible_location_indices = 1:IHT.xpoints(problem.grid);

# ╔═╡ ec444543-c871-468b-aad0-aa494586e514
md" ### refit $(@bind is_fit_on CheckBox(default = false))"

# ╔═╡ 5bc1c582-825e-4906-b796-9da5c14554fb
@bind  direct_problem_inds PlutoUI.combine() do Child
	md"""
	Select indices for direct problem 
	
	i1 = $(
		Child(Select(possible_location_indices,default=8))
	) \
	i2 = $(
		Child(Select(possible_location_indices,default=20))
	)\
	i3 = $(
		Child(Select(possible_location_indices,default=25))
	)\
	i4 = $(
		Child(Select(possible_location_indices,default=OHT.xpoints(problem.grid) - 1))
	)\
	"""
end

# ╔═╡ 6f7cd9aa-0a48-4870-9223-4b7366d06867
md" Add noise $(@bind make_noisy CheckBox(default = false)) "

# ╔═╡ 542bebdb-0cdc-4c07-9210-9c773b406bbd
md" noise amplitude $(@bind noise_amp Slider(0:100, default = 0, show_value= true))"

# ╔═╡ 99f51ecf-36e4-4a0b-aa14-0757c88c69a3
md" time region for measured data [0.0 ... $(@bind measured_tmax Slider(IHT.trange(problem.grid), default = maximum(IHT.trange(problem.grid)), show_value = true) )]"

# ╔═╡ 81314d03-95c2-4eb5-9f15-09f660b2c875
md" Lower λ limit = $(@bind lower_lam_limit confirm(Slider(0.0 : 1e-3 : 10.0, default = 0.1, show_value = true)))"

# ╔═╡ 488e4960-4420-4875-a5c4-cecba610dd1a
md" Upper λ limit = $(@bind upper_lam_limit confirm(Slider(0.1 : 1e-3 : 100.0, default = 20, show_value = true)))"

# ╔═╡ 1843f411-4706-441f-a07f-5cd77a5650b1
md" Λ(T) parameters number $(@bind basis_degree Select(1:20, default = 4))"

# ╔═╡ c1bd739f-39b0-4849-a468-cb177dab7d6a
begin 
	inds_real = unique(collect(direct_problem_inds))
	t_measured = collect(OHT.trange(problem.grid))
	flag = t_measured .< measured_tmax
	t_measured = t_measured[flag]
	measured_temperatures = transpose(problem.T[inds_real, flag])
	if make_noisy
		@. measured_temperatures *= 1 + 1e-3*noise_amp*rand()
	end
	
	p_measured = Plots.plot(grid = true,gridlinewidth=3,gridstyle = :dot,minorgrid=true, box = :on)
	for (i,c) in enumerate(eachcol(measured_temperatures))
		x_val = OHT.xvalue(problem.grid, inds_real[i])
		Plots.plot!(p_measured,t_measured, c, label = round(x_val*1e3,digits = 2), linewidth = 3)
	end
	xlabel!(p_measured, "Time, s")
	ylabel!(p_measured, "Temperature, ᵒC")
	p_measured
end

# ╔═╡ 8189333a-7ee8-44a0-9d76-a8efccd22666
md" Particle Swarm iterations number $(@bind pso_iters Select(10:10:10000 , default = 100))"

# ╔═╡ 9a602387-0383-41a3-84f2-2c31fbd93d6e
@bind refresh_plot Button("Refresh plot")

# ╔═╡ a02ce5d4-1845-410e-988c-f9f71fd129d3
@bind  lam_changing_pars PlutoUI.combine() do Child
	md"""
	### Trial value of thermal conductivity
	
	``\lambda(T) = a_0T^3 + a_1T^2 + a_2T + a_3``
	
	a0 = $(
		Child(Slider(1e-3:0.01:100,default=2,show_value = true))
	) \
	a1 = $(
		Child(Slider(1e-3:0.01:100,default=2,show_value = true))
	)\
	a2 = $(
		Child(Slider(1e-3:0.01:100,default=2,show_value = true))
	)\
	a3 = $(
		Child(Slider(1e-3:0.01:100,default=2,show_value = true))
	)\
	a4 = $(
		Child(Slider(1e-3:0.01:100,default=2,show_value = true))
	)\
	a5 = $(
		Child(Slider(1e-3:0.01:100,default=2,show_value = true))
	)\
	"""
end

# ╔═╡ 991bba10-d77b-402c-ad3b-a69c719ef58f
begin 
	time_data = t_measured
	temperatures = collect(measured_temperatures)
	initial_distribution = 20.0
	therm_locations = (i -> OHT.xvalue(problem.grid, i)).(inds_real)	
end

# ╔═╡ b7f738ad-2faa-4399-be33-71c43184a4c7
begin 
	C = IHT.OptimizableVariable(IHT.ScaledPolynomial(IHT.BernsteinSymPoly(Cp_fun.poly.coeffs), xmin = lam_T_range[1], xmax = lam_T_range[end]))
	λ_poly = IHT.ScaledPolynomial(IHT.BernsteinSymPoly(ntuple(_->1.0 , basis_degree)), xmin = lam_T_range[1], xmax = lam_T_range[end])
	λ = IHT.OptimizableVariable(λ_poly, flag = true , lb = lower_lam_limit , ub = upper_lam_limit)
	dλdT_poly =  PW.derivative(λ_poly)
	dλdT = IHT.OptimizableVariable(dλdT_poly, flag = false)

end

# ╔═╡ b7cab05e-c24f-4268-ba7b-f361c1e9076d
begin
	xpoints_number = 100
	tpoints_number = 2000
	inv_probl = IHT.SingleInverseProblem(time_data, temperatures, initial_distribution, therm_locations, C,λ, dλdT, H, xpoints_number, tpoints_number)
end

# ╔═╡ 98a565ad-cfa2-49cd-9380-d577df12cfdd
md"### Do profiling? $(@bind is_do_profile CheckBox(default = false))"

# ╔═╡ 01d72198-4ee9-436b-ada2-035821b8698d
if is_do_profile 
	@profview  discrepancy(a123,inv_probl)
end

# ╔═╡ e9c5ed77-2ea4-4c20-b1f5-88d7fe700e6c
begin 
	dir_field = IHT.temperature_field(problem)
	inv_field =IHT.temperature_field(inv_probl.direct_problem)
end;

# ╔═╡ b832d7ff-e2b8-4822-8400-7729c8ccff6f
is_fit_on || discrepancy(lam_changing_pars[1:basis_degree], inv_probl)

# ╔═╡ ba05924f-0aee-4e6c-8677-b6f10c50d5e2
if is_fit_on 
	
	(start, lb, ub)  =  IHT.fill_starting_vectors(inv_probl)
	
	fun = OptimizationFunction(IHT.discrepancy!, NoAutoDiff())
	optp = OptimizationProblem(IHT.discrepancy!, start, inv_probl, lb = lb, ub = ub , maxiters=pso_iters )
	res = solve(optp, ParticleSwarm())

end

# ╔═╡ 821e442e-85e1-446d-b468-5cb64edddae3
inv_probl

# ╔═╡ 5c4be264-6b6e-466d-a1f6-eae02dea39d1
res

# ╔═╡ a4c907f9-1817-4c84-9bd1-6c5901473103
inv_probl

# ╔═╡ 22203a92-5d3d-4416-a86d-89aff3a939cc
dλdT

# ╔═╡ ae28ced2-046a-4751-947e-9d2fb8c412f8
md"addtitional fit $(@bind is_after_fit CheckBox(default = true))"

# ╔═╡ fcd801f7-d34d-40d6-8902-8f926dc56e8b
if is_fit_on && is_after_fit
	
	new_start = res.u
	disc_before = discrepancy(new_start, inv_probl)
	
	optp2= OptimizationProblem(discrepancy, new_start, inv_probl)
	res2 = solve(optp2, NelderMead())
	disc_after = discrepancy(res2.u, inv_probl)
	disc_before - disc_after
end

# ╔═╡ aba65d0b-ae86-4005-943d-c242b335c08c
@doc Main.InverseHeatTransfer.SingleInverseProblem

# ╔═╡ e4aac0a1-d9c0-4e20-8f86-151090a9651f
#=begin 
	lam_pars
	lam_prs = Polynomials.fit(Polynomial{Float64},lam_T_range,lam_values,3)
	lam_der = Polynomials.ImmutablePolynomial( derivative(lam_prs))#;% производная теплопроводности
end;=#

# ╔═╡ 6f537a4c-bd46-465f-8846-e87a3f7d38c1
function plot_optimizable(C;kwargs...)

	return Plots.plot(range(C.p.xmax, C.p.xmin,100),C.p.(range(C.p.xmax, C.p.xmin,100)); kwargs...)
end

# ╔═╡ a302366d-27f9-4600-930d-fc0acc9456cd

begin 
	lam_changing_pars
	refresh_plot

	(Tmeas_min,Tmeas_max) = extrema(temperatures)
	(lam_minlim, lam_maxlim) = extrema( (extrema(lam_values)...,extrema(λ.(range(Tmeas_min,Tmeas_max,100)))...) )
	plot_common_args = (grid = true, gridlinewidth=3, gridstyle = :dot,minorgrid=true, box = :on, linewidth = 2)
	discr  = norm(inv_probl.residual)

	# plotting fitted values 
	p_fit_lam = plot_optimizable(λ, label = "fitted"; plot_common_args...)
	plot!(p_fit_lam, lam_T_range,lam_values, label = "real"; plot_common_args...)
	xlims!(p_fit_lam,(0.99*Tmeas_min,Tmeas_max*1.01))
	ylims!(p_fit_lam, (0.9*lam_minlim, 1.1*lam_maxlim))
	xlabel!(p_fit_lam, "Temperature, ᵒC")
	ylabel!(p_fit_lam, "Thermal conductivity, W/(m*K)")
	title!(p_fit_lam, "Fitted VS measured")

	# plotting temeprature distribution 
	
	p_distr = Plots.plot(inv_probl.Tdata_evaluated, label = "calced")
	Plots.plot!(p_distr, inv_probl.Tdata_measured, linestyle = :dash, label = "measured"; plot_common_args...)
	title!(p_distr,"discrepancy = $(discr)")
	xlabel!(p_distr, "Timestep")
	ylabel!(p_distr, "Temperature, ᵒC")
	
	
	# ploting residuals

	p_residual = Plots.plot(inv_probl.residual, label = nothing ; plot_common_args...)
	title!(p_residual, "Residuals")
	xlabel!(p_residual, "Timestep")
	ylabel!(p_residual, "ΔT, ᵒC")
end ;

# ╔═╡ 454e2649-539e-443e-bd1f-d5c2cd5b5cf6
p_fit_lam

# ╔═╡ 8c970432-b3c7-4837-9456-7744b0ffcc90
p_residual

# ╔═╡ 3736ef61-1756-4c68-95cb-1d61f57b321a
p_distr

# ╔═╡ 5ab2bb2f-fd46-42ff-bdb2-3ba08231d9b4
discr

# ╔═╡ 169ae691-dcf6-465d-ba44-5cd8c086651d
begin 
	p_fit_C = plot_optimizable(C)
	plot!(lam_T_range, cp_values)
	title!("Specific heat Cp⋅ρ")
end

# ╔═╡ 7a22e820-13df-4e56-8507-d1b089fe7f52
md" Use bench? $(@bind bench CheckBox(default = false))"

# ╔═╡ 0ed065e8-4f97-46ba-bf07-eb02fb2ae9d6
!bench || @benchmark fd_solver($problem,$solver_type)

# ╔═╡ 5e7504fd-29f1-41a3-b940-5921c471a4f5
problem

# ╔═╡ 802f203a-8f5b-4a64-a26b-4fe7f812ff56
solver_type

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

# ╔═╡ ad366a29-e978-4aa4-a3e6-a269556ca766
# Кнопка для запуска сканирования
@bind scan_trigger Button("Сканировать папку")

# ╔═╡ 5c0c8bad-d8b7-433a-a317-cfbedf4caa61
begin
	scan_trigger
	files = filter(isfile, readdir(pwd(), join=true))
	@bind selected_files confirm(PlutoUI.MultiSelect(files .= basename.(files)))


end

# ╔═╡ Cell order:
# ╠═c7400800-f152-11f0-a386-29e067f9216c
# ╠═b6133376-a592-48ee-8d63-a8413e20088d
# ╠═5fe54406-20d7-4747-bcd4-c5c0f17af47a
# ╠═e6216c49-09a8-45c5-8c3a-6d1f7281cc3e
# ╠═79ff4f54-ce5d-4f30-81ea-0c6cbdadceb0
# ╠═8a53eac7-3fa8-4d33-bf86-3cfcbdf120c5
# ╠═a4f32130-c129-426f-8092-a688e981aea8
# ╠═c806aa46-a30e-4431-b2e1-c7131ce86615
# ╠═32ac05db-67f0-4638-9842-d5ef37442efa
# ╠═1936a00f-69f3-48e4-9176-6b2186e40c2c
# ╠═6c87fcb7-55e6-44f4-8978-12e74b9cfdbf
# ╠═ea882e26-5972-443f-a7ac-3075701b90fe
# ╠═eeb3d1a0-be52-444d-aa82-89c639d5aece
# ╠═4d348c80-1035-4561-85e2-02c9bf12486a
# ╟─da943c51-3f33-4f10-bc3d-9129f76399df
# ╟─1611afbd-f13c-4532-aec9-6404b33e3f15
# ╟─d8e97e88-da66-48d7-b5cd-4c1b3d9a4c4c
# ╟─4316370d-583d-4b48-8632-d17bc34209bd
# ╟─d6145328-5f43-4148-8911-232c9fed39c4
# ╟─57b23767-880e-4375-b9b4-bd0164df4a4b
# ╟─3a88786c-cefd-4daa-bc71-1d7cf04b9a2e
# ╟─36e2c7b2-0cb0-4684-8d01-03e107349cb2
# ╟─9db7fc50-8d3e-416b-8a57-2c819d6dc426
# ╟─ec6e4678-a807-47f6-845d-dc529f38e80c
# ╟─f9f5058d-6943-4d96-9821-eafd89100fc6
# ╟─a554174d-82d6-4acd-89a9-bc21aca78e7d
# ╟─8c0b27e0-7f0b-4921-a679-ceda2ad46e82
# ╟─a4cd1c4f-f744-4ec1-bd4d-52a9e83a73ff
# ╟─e746f979-e39d-4838-8ae1-a64a1afed1f4
# ╟─ed45f038-ae2b-4888-8547-ab9638170842
# ╟─ddc9ba47-c4f2-4c77-8a9c-c598b7eb322b
# ╟─3acc0be7-dabc-4577-917d-26530cf192bd
# ╟─38976547-9399-49f1-b83e-0833327ea0ed
# ╟─8f6d3be5-ffb5-41d5-8e50-f4353193c53c
# ╟─0eae0fd4-8b5d-448a-b935-e06d469f080b
# ╟─9fbe7c0d-90b7-4348-8b22-d8058452a02b
# ╟─84bfb85d-def3-4a08-be60-01de2d68be36
# ╟─8992ffda-fcbb-4994-9acb-1080056903b7
# ╟─331ee2df-931b-4267-8795-1fe402c7fdaa
# ╟─86a00c2e-4896-417d-93d0-bc69b5823957
# ╟─95508f3b-73b5-4cc5-b6bd-bb4d949bf00c
# ╟─1182c47f-3cfb-4aed-a597-c79a0a786d3c
# ╟─ec444543-c871-468b-aad0-aa494586e514
# ╟─5bc1c582-825e-4906-b796-9da5c14554fb
# ╟─6f7cd9aa-0a48-4870-9223-4b7366d06867
# ╟─542bebdb-0cdc-4c07-9210-9c773b406bbd
# ╟─99f51ecf-36e4-4a0b-aa14-0757c88c69a3
# ╟─81314d03-95c2-4eb5-9f15-09f660b2c875
# ╟─488e4960-4420-4875-a5c4-cecba610dd1a
# ╟─1843f411-4706-441f-a07f-5cd77a5650b1
# ╟─c1bd739f-39b0-4849-a468-cb177dab7d6a
# ╟─8189333a-7ee8-44a0-9d76-a8efccd22666
# ╠═454e2649-539e-443e-bd1f-d5c2cd5b5cf6
# ╟─9a602387-0383-41a3-84f2-2c31fbd93d6e
# ╟─a02ce5d4-1845-410e-988c-f9f71fd129d3
# ╠═8c970432-b3c7-4837-9456-7744b0ffcc90
# ╠═3736ef61-1756-4c68-95cb-1d61f57b321a
# ╟─991bba10-d77b-402c-ad3b-a69c719ef58f
# ╠═b7f738ad-2faa-4399-be33-71c43184a4c7
# ╠═b7cab05e-c24f-4268-ba7b-f361c1e9076d
# ╠═98a565ad-cfa2-49cd-9380-d577df12cfdd
# ╠═01d72198-4ee9-436b-ada2-035821b8698d
# ╠═e9c5ed77-2ea4-4c20-b1f5-88d7fe700e6c
# ╠═5ab2bb2f-fd46-42ff-bdb2-3ba08231d9b4
# ╠═b832d7ff-e2b8-4822-8400-7729c8ccff6f
# ╠═a302366d-27f9-4600-930d-fc0acc9456cd
# ╠═ba05924f-0aee-4e6c-8677-b6f10c50d5e2
# ╠═821e442e-85e1-446d-b468-5cb64edddae3
# ╠═5c4be264-6b6e-466d-a1f6-eae02dea39d1
# ╠═a4c907f9-1817-4c84-9bd1-6c5901473103
# ╠═22203a92-5d3d-4416-a86d-89aff3a939cc
# ╠═ae28ced2-046a-4751-947e-9d2fb8c412f8
# ╠═fcd801f7-d34d-40d6-8902-8f926dc56e8b
# ╟─169ae691-dcf6-465d-ba44-5cd8c086651d
# ╠═aba65d0b-ae86-4005-943d-c242b335c08c
# ╟─e4aac0a1-d9c0-4e20-8f86-151090a9651f
# ╠═6f537a4c-bd46-465f-8846-e87a3f7d38c1
# ╟─7a22e820-13df-4e56-8507-d1b089fe7f52
# ╟─0ed065e8-4f97-46ba-bf07-eb02fb2ae9d6
# ╟─5e7504fd-29f1-41a3-b940-5921c471a4f5
# ╟─802f203a-8f5b-4a64-a26b-4fe7f812ff56
# ╟─7017c716-8682-4a74-9bfe-7dc09879db16
# ╟─e525b2ec-f4a7-4052-af2b-7a9b4fa4f520
# ╟─ad366a29-e978-4aa4-a3e6-a269556ca766
# ╠═5c0c8bad-d8b7-433a-a317-cfbedf4caa61
