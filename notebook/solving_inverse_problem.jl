### A Pluto.jl notebook ###
# v1.0.1

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

# ╔═╡ a17fe1fe-5542-454b-b45e-942ac52b6f1a
begin 
	import Pkg
	Pkg.activate(@__DIR__)
	Pkg.instantiate()
	using Interpolations  
	using StaticArrays , Plots  , Revise , PrettyTables
	using Optimization , OptimizationOptimJL , OptimizationNLopt , OptimizationMetaheuristics
	using PlutoUI ,  OrderedCollections
	using HypertextLiteral
	using Revise
	import InverseHeatTransfer
	using Dates
	using Tables
	using CSV
	using FiniteDiff
	using LinearAlgebra
	using StatsBase, BenchmarkTools
	using Distributions
end

# ╔═╡ 5807712b-5d26-49c8-ab65-dac167ebad7b
begin 

	data_selection_path_ref = Ref("") # data selection saving path
	data_selection_name_ref = Ref("") # data selection saving name
	
	default_data_fodler = joinpath(@__DIR__, "..","src","data_utils","property_inversion_ansys_new")
	source_path = joinpath(@__DIR__,"..","src")
	includet(joinpath(@__DIR__ , "CustomPlutoFunctions.jl"))
	
	IHT = InverseHeatTransfer
	PW = IHT.ScaledPolynomials
	OHT = IHT.OneDHeatTransfer
	DC = IHT.DataConnector
	
	plot_common_args = (grid = true, gridlinewidth=3, gridstyle = :dot,minorgrid=true, box = :on, linewidth = 3)
end;

# ╔═╡ 2bfb4e52-6248-4832-aca1-98ba58959bff
PF = Main.CustomPlutoFunctions;

# ╔═╡ 3b1c3b0a-558e-4987-bf16-072963e455cf
html"""
<style>
.plutoui-rangeslider{
		width : 100%;
}
.plutoui-slider{
		width : 100%;
}
</style>
"""

# ╔═╡ 2db7cd7e-3d1a-4533-b0de-27188300565b
const OPTIMIZERS = [
OptimizationMetaheuristics.PSO => "Particle swarm", OptimizationNLopt.NLopt.LN_COBYLA => "COBYLA", 
OptimizationNLopt.NLopt.LN_BOBYQA => "BOBYQA" , 
OptimizationOptimJL.LBFGS => "LBFGS", 
OptimizationNLopt.NLopt.LD_LBFGS =>"NL-LBFGS"
];

# ╔═╡ db671921-13dc-497b-81e5-dcb4da0695f9
md""" ## Data selection file loading """

# ╔═╡ 450fb200-eec6-4e96-9ebd-81453c015830
md" ##### Load data from : $(@bind input_data_type Select([:hdf5_data_selector] , default = :winpos))"

# ╔═╡ 6e062bd9-d20c-4e1d-b772-328bec8859ea
begin 
	md""" 
	
	#### working folder $(@bind data_selection_folder TextField(90, default = realpath(default_data_fodler))) 
	
	"""
end

# ╔═╡ a7abe643-2553-450d-ac81-d4a690a1c2ff
if isdir(data_selection_folder) 
	
	all_hdf5_files = [d for d in readdir(data_selection_folder) if contains(d , ".hdf5") ]
	if isempty(all_hdf5_files)
		md" **There is no hdf5 files in the folder**"
	else
		md" **Select experiment file** $(@bind data_selection_name Select(all_hdf5_files))"
	end
else
	md" Incorrect folder $(data_selection_folder)"
end

# ╔═╡ 41bc1a0a-73c8-430d-a1d3-4eb98487c815
@bind reload_trigger Button("Reload")

# ╔═╡ 34b63c47-678d-4a5f-aed3-1896b778e117
try
	fful_name = joinpath(data_selection_folder , @isdefined(data_selection_name) ? 	data_selection_name : "" )
		Text(DC.read_measurements_specification_string(fful_name ,  specification_node_key = "x_comment"))
	catch er 
		md" There is no file  $(er)"

end

# ╔═╡ 444236e5-e010-4b8b-8709-31c0307dc5d8
begin 
	reload_trigger
	data_selection_path_ref[] = data_selection_folder
	@isdefined(data_selection_name) && (data_selection_name_ref[] = data_selection_name)
end;

# ╔═╡ 81c2f93b-05b1-4eb0-9919-4ef76ecad233
begin 
	reload_trigger
	is_data_loaded = false
	all_data = nothing
	data_selection_fullfile = joinpath(data_selection_path_ref[] , data_selection_name_ref[])
	try 
	
		if input_data_type == :hdf5_data_selector
			global all_data = DC.WinPos.load_from_hdf5(data_selection_fullfile , DC.DataSelectorsGroup)
			global is_data_loaded = true
		else
			global is_data_loaded = false
			error("Incorrect data type")
		end
	catch er
		
		global is_data_loaded = false
		md" Not loaded $(er)"
	end
	is_data_loaded = !isnothing(all_data)
end;

# ╔═╡ 074c47c6-a1ae-4ded-8e64-b40393d7ba4a
!is_data_loaded && md" there is no such file $(data_selection_fullfile)"

# ╔═╡ 35b5c27e-921f-4fe7-90bc-3b327d9150fe
if is_data_loaded 
	
	
	table_data = Matrix{Any}(undef , (length(all_data.d) , 6))
	for (i , (k , d)) in enumerate(all_data.d)
		table_data[i, 1]  = "P$(i)"
		table_data[i, 2]  = k
		table_data[i, 3] = 1e3*DC.thickness(d)
		table_data[i, 4] = [ Pair(v,k) for (k,v) in zip(1e3*DC.sensors_locations(d) , DC.selected_names(d))]
		table_data[i, 5] = DC.tmin(d)
		table_data[i, 6] = DC.tmax(d)
	end

	input_table = pretty_table(HTML , table_data , column_labels = ["probl", "name"," h" , "Locs","tmin", "tmax"] , top_left_string ="Sample properties")

end

# ╔═╡ 5be15388-cea2-4884-b8de-bff5be64e506
if is_data_loaded	
	PN = length(all_data.d)
	projects_names = ["P$(i)" for i in 1:PN]
	
	reload_trigger
	
	raw_data_plot = Plots.plot(;plot_common_args...)
	for (p_n , (k , d_i))  in enumerate(DC.selected_data_cutted_with_keys(all_data))
		isempty(d_i.data) && continue
		_t = d_i.data[:,1]
		_names = d_i.names
		CN = size(d_i.data, 2)
		for (i, c) in enumerate(eachcol(d_i.data)[2:CN])
			Plots.plot!(raw_data_plot , _t , c , label ="P$(p_n) : $(_names[i + 1])" , linestyle = :auto ; plot_common_args...)
		end
	end
	raw_data_plot
	xlabel!(raw_data_plot , "Time, s")
	ylabel!(raw_data_plot , "Temperature, °C")
else
	PN = 0
	projects_names = [""]
end

# ╔═╡ 6c4cb363-3bda-4dfa-8499-8201a013895d
is_data_loaded && @bind show_data_table Select(collect(keys(all_data)))

# ╔═╡ b31a7c52-7c4b-480f-84d8-fcb1c71dca1b
if is_data_loaded
	reload_trigger
	
	selected_plot = Plots.plot(;plot_common_args...)
	_data_combined = DC.combine_selected_data(all_data[show_data_table])
	_t = _data_combined.time_data
	for (i,c) in enumerate(eachcol(_data_combined.temperatures))
		Plots.plot!(selected_plot , _t , c ; label=_data_combined.selected_names[i]  , plot_common_args...)
	end
	title!(selected_plot , show_data_table)
	xlabel!(selected_plot , "Time , s")
	ylabel!(selected_plot , "Temperature , oC")
	selected_plot
end

# ╔═╡ c69d07a3-ba0a-48a2-b296-1944b8cd322c
is_data_loaded && @htl("""
<div style="max-height: 300px; overflow-y: auto; border: 1px solid #ccc;">
    $(pretty_table(HTML , hcat(_data_combined.time_data, _data_combined.temperatures) , column_labels = ["t" , _data_combined.selected_names...] , top_left_string ="Temeratures for $(show_data_table)"))
</div>
""")

# ╔═╡ 62069546-44fb-4a77-986d-f03624719e29
md" ## Physical properties"

# ╔═╡ a46970c8-7c91-4795-83a9-41c56a7ca399
md""" ### Compare to reference $(@bind show_passport CheckBox(default = true))"""

# ╔═╡ 74707b89-d91c-468d-9b4a-a06dc99f69c0
md"### Heat capacity setup"

# ╔═╡ 4f3e1899-541d-4809-810c-f2e6b7ca1ad1
md"#### Optimize C : $(@bind is_optimize_c CheckBox(default = false))"

# ╔═╡ 4fdbfb71-ec8d-4d30-b7e7-f289f773fadd
md" **C(T) parameters number** $(@bind c_basis_degree Select(1:100, default = 4))"

# ╔═╡ eb755603-1eda-4182-a026-e9b163e78ae3
md" **Use table for constraints** $(@bind use_cp_table_for_constraints CheckBox(false))"

# ╔═╡ 97618f0c-52ae-4772-b5b7-5512ea44af09
!use_cp_table_for_constraints  && md" **Upper Cp limit** = $(@bind upper_cp_limit confirm(Slider(100.0 : 1e-1 : 10000.0, default = 1300.0, show_value = true)))"

# ╔═╡ baffb68a-fb52-4bac-b484-4a215108aaed
if !use_cp_table_for_constraints 
	md" **Lower Cp limit** = $(@bind lower_cp_limit confirm(Slider(100.0 : 1e-1 : 10000.0, default = 750.0, show_value = true)))"
end

# ╔═╡ 8345302e-0b9d-4208-9a66-3d9f32903b39
md""" Show Cp confidence bounds $(
		@bind is_show_cp_conf CheckBox(true)
	)
"""

# ╔═╡ 0ff89750-2a97-436b-b759-7da1354b2c6f
md" ## Thermal conductivity setup"

# ╔═╡ 4ca95124-8a7a-4e4f-9e65-ff7b2adf35a5
md" #### Optimize λ ? $( @bind is_optimize_lambda CheckBox(default = false))"

# ╔═╡ fa72774e-040a-4bc3-a759-5eb68c243fb4
md" Λ(T) parameters number $(@bind lam_basis_degree Select(1:100, default = 4))"

# ╔═╡ 88e8d37d-e4e4-486d-921e-03a74fbf00f2
md" **Use table for constraints** $(@bind use_lam_table_for_constraints CheckBox(false))"

# ╔═╡ c8dc4f9d-a549-4dc2-82bd-38ffe949ea55
!use_lam_table_for_constraints && md"**Lower λ limit =** $(@bind lower_lam_limit confirm(Slider(0.0 : 1e-3 : 10.0, default = 0.1, show_value = true)))"

# ╔═╡ bfa23359-8bac-4db0-bac1-0885ebe8ec4b
!use_lam_table_for_constraints && md" **Upper λ limit =** $(@bind upper_lam_limit confirm(Slider(0.1 : 1e-3 : 30.0, default = 20, show_value = true)))"

# ╔═╡ 68ed85b2-7308-4ea7-b696-0f1951219592
@bind lam_y_scale_region PlutoUI.combine() do Child 
md"""
	**Figure y-range**
	
	``\lambda_{min} ``= $(
		Child(NumberField(0:1e-2:50,default=0.0))
	) \
	``\lambda_{max} `` = $(
		Child(NumberField(0:1e-3:50,default=30.0))
	)\

	show confidence bounds $(
		@bind is_show_lambda_conf CheckBox(true)
	)
	"""
end

# ╔═╡ 9b35c13b-24ba-4c43-a621-a3f8ba45fe4a
md"""
--------------------
### Inverse problem setup

"""

# ╔═╡ fc5c1209-26ec-41d7-9238-e57d18330de1
@bind refit Button("Refit!!")

# ╔═╡ b1f00c55-5ce9-4a0f-a548-a8a9041d02fc
if is_data_loaded 
	is_data_ready = true
	try 
		refit
		for (k , d_i) in all_data.d
		 	DC.combine_selected_data(d_i)
		end
		global is_data_ready = true
	catch 
		global is_data_ready = false
	end
	is_data_ready
end

# ╔═╡ 2b3a1a65-8d3c-424e-a6cf-0e96646795f4
md" ### Refit $(@bind is_fit_on CheckBox(default = false))"

# ╔═╡ 3432fabd-911b-4370-af86-c396cbc7bbab
md" #### Individual settings for each problem $(@bind is_individual_settings CheckBox(false))"

# ╔═╡ 33473b7a-e22f-4333-a2f2-374778c0d603
md"""



Regularization multiplier α: $(@bind reg_multiplier confirm(NumberField(0.0 : 1e-3 : 1000 , default = 1e-3)))
"""

# ╔═╡ aad5d954-fb91-4a7d-a09e-02105ef3d0f9
md"""
------------
### Optimization results
"""

# ╔═╡ 5e5b27d1-fdc8-451d-a4f0-93f33adbcf76
md" Display residuals $(@bind residuals_type Select([:weighted , :unweighted , :both] , default = :both))"

# ╔═╡ 5d4b49c0-1d0f-41b5-b359-30c033fc8544
input_table

# ╔═╡ 2a1ca349-3732-443e-bc48-5be611a5d91f
md"""
### Select optimizer $(@bind optimizer Select(OPTIMIZERS, default = OptimizationNLopt.NLopt.LN_BOBYQA)) 
"""

# ╔═╡ d02ec3fd-1b0d-4bc3-92e9-adedc8bf2c8c
md" ##### Global optimizer iterations number $(@bind pso_iters Select(10:10:10000 , default = 200))"

# ╔═╡ 0a0324df-430f-4ac0-857b-4da6a7dca138
md"addtitional fit $(@bind is_after_fit CheckBox(default = false))"

# ╔═╡ c8cd1797-b631-4fdf-a51f-67e22c86c55a
md"""
	### L-curve analysis $(@bind use_l_curve CheckBox(default = false))

	α points number $(@bind l_curve_points_number Select(3:100 , default = 20) )
	"""

# ╔═╡ aea6fc8f-5ee2-40a3-aebd-942c45eec0d6
md" Use logarithmic scale $(@bind is_l_curve_logscale CheckBox( default = true))"

# ╔═╡ a8f5cc07-b2bd-445c-8797-ffe78a25641b
begin 
	if is_l_curve_logscale
		l_culve_alphas_default_range = range(-4.0 , 2 , 1000)
	else
		l_culve_alphas_default_range = range(1e-6 ,10 , 1000)
	end
end

# ╔═╡ eb70f6d9-301f-4e0f-91a2-0e911b1a6dd9
md"""

 log10(α) search range: $(@bind l_curve_alpha confirm(RangeSlider(l_culve_alphas_default_range)) )


"""

# ╔═╡ 10f7c083-a697-4b7e-87ea-de2791ed1a30
md" α = $(10 .^extrema(l_curve_alpha)) " 

# ╔═╡ fc79d398-03d0-4f75-a777-41e2e27fee5e
if is_data_loaded 
	lam_T_range = DC.default_temperature_range(all_data)
	(Tmin,Tmax) = (lam_T_range[1],lam_T_range[2])
	Tplot = range(Tmin, Tmax, 100)
end;

# ╔═╡ ffeafacd-2b98-48c0-b840-297fb51f5e54
input_table

# ╔═╡ ea412a80-0bf7-4ac6-9b90-8530a4c26008
md"### Save data to file ? $(@bind is_write_file CheckBox(false))"

# ╔═╡ 23902c7d-5973-4868-8488-e0c7634573c4
@bind default_save_file TextField(60,default=first(eachsplit(data_selection_name , "."))*"_solution.hdf5")

# ╔═╡ 6800ae42-c5b1-4d1e-84fb-a03015bf138f
@bind save_folder TextField(90, default = realpath(default_data_fodler))

# ╔═╡ c4053281-7cd4-4590-b004-b4ca43771094
md"""### Click to resave file $(@bind resave_file Button("resave"))"""

# ╔═╡ 9cd8c4c8-7d11-4fb0-92d5-439702aa9496
md" ### OPTIMIZATION "

# ╔═╡ 70e2c8f2-d896-4773-8280-d391d9975307
function extract_diag(m::AbstractMatrix)
	N = size(m , 1)
	@assert N==size(m, 2) "Matrix must be square"
	return [m[i , i] for i in 1:N]
end

# ╔═╡ bce5f7ae-49c5-4520-a0f6-6215e5078674
includet(joinpath(source_path, "problem_ensemble_functions.jl"))

# ╔═╡ d71fd6d1-167d-40fb-a253-b0982e19c0d0
@bind Ttest Slider(range(lam_T_range[1],lam_T_range[2],100), show_value = true)

# ╔═╡ 58e0b175-a114-46be-9792-35954faea5ec
md" ### MATERIAL PROPERITES "

# ╔═╡ efa9120f-5a45-41a0-9132-dc26f967fec3
begin 
	Rho_Dict = Dict(
	
	"RBSN" => 2700.0, 
	"NIASIT" => 2000.0,
	"OTM357" => 2530.0,
	"HAFS" => 1600.0,
	"BAFS" => 1800.0
	
);
	md" ### ρ"
end
	

# ╔═╡ fb2937d0-738b-4329-aef5-f3af1d17497a
begin 
	Cp_Dict = Dict(
	"RBSN" =>(;
		x = Float64.([ 20.0 ,100, 500, 800, 1100, 1500, 2000, 2500]),
		y = Float64.([690, 810, 1160, 1240, 1250, 1300, 1350, 1450])
			 ), 
	"NIASIT" => (;
		x =  Float64.([20, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100]),
		y = 1000.0 * [0.74, 0.85, 0.96, 1.04, 1.09, 1.14, 1.17, 1.2, 1.23, 1.25, 1.27, 1.29]
				),
	"OTM357" =>(;
		x = Float64.([50,100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100]),
		y = 1000.0 * [0.81, 0.92, 1.03, 1.09, 1.14, 1.17, 1.2, 1.22, 1.24, 1.26, 1.27, 1.27]
		),
	"BAFS" => (;
		x = Float64.([20, 100,200,300,400,500,600,700,800]),
   		y = Float64.([730,800,920,1000,1050,1080,1090,1090,1090])
			 )
);
	md" ### Cₚ"
end

# ╔═╡ 8f65e989-1abe-4689-9c3f-fdb6cabd7eae
md" #### Select material name : $(@bind material_name Select(collect(keys(Cp_Dict))))"

# ╔═╡ 7d37019a-34bf-43ca-aa81-8b6c29191473
md"""

### Density, kg/m³ $(@bind density confirm(NumberField(200.0:0.1:4000, default = Rho_Dict[material_name])))

"""

# ╔═╡ 5843fcfb-4e0e-480d-b263-e3f3ad6a7ac3
if density != Rho_Dict[material_name]
	md" passport density is $(Rho_Dict[material_name])"
end

# ╔═╡ 023dc25f-6cf8-4802-83b4-77d1827cd2a7
begin 
	L_Dict =Dict( "RBSN" =>(;
						x = Float64.([ 25, 100, 300, 500, 700, 900, 1100, 1300, 1500, 1600, 1700, 1800, 2000, 2500]),
						y = Float64.([18, 16, 14,  12,  11,  9,    8,    7,     7,   6,     7,    11, 12, 14])
					   ),
									   
"NIASIT" =>(;
	x = Float64.([20, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100]),
	y = [0.83, 0.88, 0.92, 0.95, 0.98, 1.01,1.05, 1.09,1.13,1.17, 1.2, 1.23]
	),
"OTM357" => (;
x = Float64.( [50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100]),
y = [1.72, 1.73,1.74, 1.75, 1.76, 1.77, 1.78,1.79,1.8, 1.82, 1.83, 1.84]
			),
"BAFS" => (;
	x = Float64.([20, 100,200,300,400,500,600,700,800]),
	y = [0.54,0.57,0.61,0.65,0.65,0.66,0.65,0.65,0.63]
	)
);
	md" ### λ"
end

# ╔═╡ ab211f5e-9cf8-438c-a3cc-4fa18960b67b
md" ### INTERNAL FUNCTIONS"

# ╔═╡ 2bf91a7e-0da8-48e8-a779-962be2e7c03b
function plot_optimizable(C; xmin,xmax , kwargs...)
	x = range(xmin , xmax ,100)
	return (Plots.plot(x, C.(x); kwargs...), x)
end

# ╔═╡ fd51a6c8-6569-4bfd-86ac-883c648fe6d9
function filter_couple_breakage!(T)
	N = size(T,1)
	for (i,c) in enumerate(eachcol(T))
		inds = findall(c .<= -196.0)
		isempty(inds) && continue
		for j in inds
			left = findfirst(Base.Fix2(>=, -196.0), view(c,j:-1:1))
			right = findfirst(Base.Fix2(>=, -196.0), view(c,j:N))
			if (isnothing(left) || isnothing(right))
				if isnothing(left) && isnothing(right) 
					c[j] = 0.0
				elseif isnothing(left)
					c[j] = c[j + right - 1]
				else
					c[j] = c[j - left + 1]
				end
				continue
			end
			tleft = c[j - left + 1]
			tright = c[j + right - 1]
			c[j] = 0.5*(tleft + tright)
		end
	end
	return T
end

# ╔═╡ dbd6f9d9-4b44-4238-86ca-7cd361545a56
function prepare_constraints_gui(el_number ,  lam_T_range ; 
					whole_range::NTuple{2 , Float64} = (600.0 , 2000.0) ,   	def_range::NTuple{2 , Float64} = (900.0 , 1300.0)  , npoints::Int = 1000 ,  
					title= "no title")
	
	_bern_max = IHT.ScaledPolynomials.scale_ξ_to_x(  		IHT.ScaledPolynomials.bern_max_locations(IHT.BernsteinSymPoly{el_number , Float64}()),lam_T_range[1] , lam_T_range[2])

	_names = ntuple(el_number) do i 
		"T=$(string(round(Int , _bern_max[i])))"
	end
		
	_default_values = ntuple(el_number) do _ 
		range(whole_range... , npoints)
	end
	_defaults = ntuple(el_number) do i 
		d = _default_values[i]
		_i =( findlast(t->t<=first(def_range), d) , findfirst(t->t>=last(def_range) , d))
		d[first(_i) : last(_i)]
	end
	PF.multi_values_table(PlutoUI.RangeSlider , _names ; title = title, default_values = _default_values  , defaults =  _defaults)
end

# ╔═╡ c8cda804-9a2a-40c9-aa91-7a48500e85ff
if  use_cp_table_for_constraints 
	@bind cp_limit_table confirm(prepare_constraints_gui(c_basis_degree , lam_T_range, title = "Cp constraints"))
end

# ╔═╡ 4330cdcd-24fc-459b-8d29-a935c6a7c347
if  use_lam_table_for_constraints 
	@bind lam_limit_table confirm(prepare_constraints_gui(lam_basis_degree  , lam_T_range; title = "λ constraints" , whole_range=(0.01,50.0) , def_range=(2.0 , 20.0) , npoints = 10000))
end

# ╔═╡ ce3dc022-0f15-41cb-8cd2-2c33b726c482
(lam_lower_bounds , lam_upper_bounds) = if use_lam_table_for_constraints 
	(l_l = ntuple(lam_basis_degree) do i 
		minimum(lam_limit_table[i])
	end,
	u_l = ntuple(lam_basis_degree) do i 
		maximum(lam_limit_table[i])
	end)	
else
	(lower_lam_limit , upper_lam_limit) 
end;

# ╔═╡ c8861fa2-9cb3-4349-9ece-eba285c24eab
if is_data_ready 
	λ_poly = IHT.ScaledPolynomial(IHT.BernsteinSymPoly(ntuple(_->1.0 , lam_basis_degree)), xmin = lam_T_range[1], xmax =lam_T_range[2])
	
	if is_optimize_lambda
		
		λ = IHT.OptimizableVariable(λ_poly, 
									flag = is_optimize_lambda , 
									lb = lam_lower_bounds , 
									ub = lam_upper_bounds )
		
		dλdT_poly =  IHT.ScaledPolynomials.derivative(λ_poly)
		dλdT = IHT.OptimizableVariable(dλdT_poly, flag = false)
		
	else
		l_x = L_Dict[material_name].x 
		l_y = L_Dict[material_name].y
		fl = @. (l_x <= lam_T_range[2]) & (l_x >= lam_T_range[1])
		l_x = l_x[fl]
		l_y = l_y[fl]
		#Cp = linear_interpolation(cp_x ,density.*cp_y , extrapolation_bc=Line())
		λ = PW.polyfit!(λ_poly, l_x, l_y)
		dλdT = PW.derivative(λ)
	end
end;

# ╔═╡ 206e1e8f-16f3-4144-8401-c2cbf40c0125

function bounds_convert(b::T , c::Union{T , Nothing} = nothing) where T <: Number
		_c = isnothing(c) ? one(T) : c
		return b * _c
end

# ╔═╡ 4b8ddfc0-ed4f-4b66-b752-fa075d348608
function bounds_convert(b::Union{Tuple , NamedTuple , AbstractArray} ,  c::T = 1.0) where T<: Number 
	return ntuple(length(b)) do i 
		c * b[i]
	end
end

# ╔═╡ d051f5e5-f836-4043-9326-639b40acaf87
(cp_lower_bounds , cp_upper_bounds) = if use_cp_table_for_constraints 
	l_cp = ntuple(c_basis_degree) do i 
		minimum(cp_limit_table[i])
	end
	u_cp = ntuple(c_basis_degree) do i 
		maximum(cp_limit_table[i])
	end	
	(
		bounds_convert(l_cp , density) , 	
		bounds_convert(u_cp , density)
	)
else
	(ntuple((_)->density*lower_cp_limit , c_basis_degree) , ntuple((_) -> density*upper_cp_limit  , c_basis_degree)) 	
end;

# ╔═╡ 16337bb4-438e-4b86-bdc5-b88ef210a960
begin 
	C_poly = IHT.ScaledPolynomial(IHT.BernsteinSymPoly(ntuple(_->density*1000.0 , c_basis_degree)), xmin = lam_T_range[1], xmax =lam_T_range[2])
	c_x = Cp_Dict[material_name].x 
	c_y = density.*Cp_Dict[material_name].y
	flc = @. (c_x <= lam_T_range[2]) & (c_x >= lam_T_range[1])
	c_x = c_x[flc]
	c_y = c_y[flc]
	if is_optimize_c
		C = IHT.OptimizableVariable(C_poly, flag = is_optimize_c , lb = cp_lower_bounds , ub = cp_upper_bounds )
	else

		#C = PW.polyfit!(C_poly, c_x,c_y)#linear_interpolation(c_x ,c_y , extrapolation_bc=Line())
		C =linear_interpolation(c_x ,c_y , extrapolation_bc=Line())
	end
end;

# ╔═╡ 04ad22b9-e474-41ed-bc87-7dbf9a2b1dcc
function eval_conf_bounds(st , ov; α = 2.0)
	@assert length(st.std) == IHT.ScaledPolynomials.parnumber(ov) "Inappropriate size"
	@warn all(i->i > 0 , st.std) "All elements of std must be greater then zero"
	
	c = IHT.ScaledPolynomials.coeffs(ov)
	v_up = deepcopy(ov)
	v_low = deepcopy(ov)
	IHT.ScaledPolynomials.refill!(v_up , @. c + α*sqrt(abs(st.std)))
	IHT.ScaledPolynomials.refill!(v_low , @. c - α*sqrt(abs(st.std)))
	return (;up = v_up , low = v_low)
end

# ╔═╡ 3b9e739c-9519-4efa-b2da-cac5451d55d3
function covariance(p , u)
	probs = deepcopy(p)
	IHT.discrepancy!(u , probs)
	f(x) = IHT.discrepancy!(x , deepcopy(probs))
	N = sum([length(p.residual) for p in probs.problems])
	σ = sum(IHT.evaluate_loss , p.problems)/N#sum([sum(t->t^2 , p.residual) for p in probs.problems])/(N^2)
	Σ  = inv(FiniteDiff.finite_difference_hessian(f , u ))
	#@show Σ
	Cov = Σ * σ
	s = extract_diag(Cov)
	return (;std = s , Cov=Cov , Σ = Σ , σ=σ , N=N)
end

# ╔═╡ f799c1f9-ebc2-4fdb-839e-ee721ab9b8f5
function descriptive_stats_table(p)
	ips = IHT.IPstats(p)
	fn = fieldnames(IHT.IPstats)
	N  = length(fn)
	tbl = Matrix{Any}(undef, (N , 2) )
	for i in 1:N 
		f_cur = fn[i]
		tbl[i , 1] = String(f_cur)
		tbl[i , 2] = getfield(ips , f_cur)
	end

	return pretty_table( 
		HTML , 
		tbl ,
		column_labels = ["name", "value"] ,
		title ="Simple descriptive statistics")
end

# ╔═╡ 544c4618-5b5e-4b25-a916-a43fc2e05913
function sensitivity_stats_table(probs , u)
	out = IHT.sensitivity_analysis_statistics(probs , u)
	tbl = Matrix{Any}(undef, (length(out) , 2))
	for (i , f_i) in enumerate(pairs(out)) 
		tbl[i , 1] = "$(first(f_i))-optimality" 
		tbl[i , 2] = last(f_i) 
	end
	
	return pretty_table( 
		HTML , 
		tbl ,
		column_labels = ["name"  , "value"] ,
		title ="Optimality criteria")
end

# ╔═╡ e83d394a-f06f-48de-a664-aad82e4361fc
function regression_stats_table(stats , u; α = 1.96)
	std_d = stats.std
	N = length(std_d)
	@assert N == length(u) "Incorrect u size"
	tbl = Matrix{Any}(undef, (N , 4) )
	for i in 1:N 
		tbl[i , 1] = "β$(i)"
		tbl[i , 2] = u[i]
		tbl[i , 3] = std_d[i]
		tbl[i , 4] = α * std_d[i]
	end

	return pretty_table( 
		HTML , 
		tbl ,
		column_labels = ["name","value" ,  "std", "d"] ,
		title ="Optimization variables ")

end

# ╔═╡ baaf2f95-9529-4829-a8fa-0724fd21b78e
function autocorrelation_stats_table(autocor_stats; is_unweighted::Bool = false)
	N = length(autocor_stats) # number of problems 
	a1 = autocor_stats[1][1]
	fn = filter(f->isa(getfield(a1 , f) , Number) , fieldnames(typeof(a1))) # tests_names 
	M  = length(fn)
	tbl_vect = Any[]
	for (i , a_i) in enumerate(autocor_stats)
		# a_i - problem stats (vector, each element - corresponds to )
 	   	TPnumber = length(a_i) # length of  	
		tbl_i = Matrix{Any}(undef, (TPnumber ,  M + 1) ) # i'th thermocouple table
		for tp_i in 1 : TPnumber #over couples
			a_i_t = a_i[tp_i] #
			row = Any[]
			push!(row, "P$(i) : T$(tp_i)")
			for (j , f_cur) in enumerate(fn)
				  push!(row , getfield(a_i_t, f_cur))
			end
			tbl_i[tp_i , :] = row
		end
		push!(tbl_vect , tbl_i)
	end
	col_names = Vector{String}(undef,  M + 1)
	col_names[1] = "names"
	for i in 1 : M 
		col_names[i + 1] = "$(fn[i]) "
	end	
	tbb = vcat(tbl_vect...)
	return pretty_table(HTML,tbb  , column_labels=col_names ,  title = "Autocor. tests for $(is_unweighted ? "unweighted res." : "weighted res.")") 
end

# ╔═╡ b32b35a1-f31c-428f-a9ab-d640f2be496e
function result_table(res)
	isnothing(res) && return nothing
	fn = fieldnames(typeof(res.stats))
	tbl = Matrix{Any}(undef , (length(fn) + 3, 2))
	tbl[1 , 1] = "u"
	tbl[1 , 2] = res.u
	tbl[2 , 1] = "objective"
	tbl[2 , 2] = res.objective
	for (i , f_i) in enumerate(fn) 
		tbl[i + 2 , 1] = f_i 
		tbl[i + 2, 2] = getfield(res.stats , f_i)
	end
	tbl[end , 1] = "Original"
	tbl[end , 2] = res.original
	pretty_table(HTML , tbl ,  column_labels=["name" , "value"] , title="Optimizer statistics")
end

# ╔═╡ 9c4936fd-25d6-4760-85f5-38aba54f1800
function rep_tup(a , N::Int)
	return ntuple(N) do _
		a
	end
end

# ╔═╡ 3281dd3c-0453-4098-a6ed-4d771717033b
md""" 

##### Direct problem number of coordinate steps 
$(@bind dp_xpoints if !is_individual_settings
	 Select(10:2000, default = 140)
else
	confirm(PF.multi_values_table(PlutoUI.Select , projects_names , default_values =rep_tup(10:2000 , PN) , defaults = rep_tup(140 , PN) , title = ""))
end)
"""

# ╔═╡ 63f66c3f-1489-4da1-963b-47812ca245ae
md""" 

##### Direct problem number of time steps 
$(@bind dp_tpoints if !is_individual_settings
	 Select(10:10:5000, default = 2000)
else
	confirm(PF.multi_values_table(PlutoUI.Select , projects_names , default_values =rep_tup(10:2000 , PN) , defaults = rep_tup(2000 , PN) , title = ""))
end)
"""

# ╔═╡ 107318d8-051f-4967-8dc0-49d836cd9469
md""" 

##### Covariance type
$(@bind covariance_type if !is_individual_settings
	 Select(IHT.ALL_COVARIANCE_TYPES)
else
	confirm(PF.multi_values_table(PlutoUI.Select , projects_names , default_values =rep_tup(IHT.ALL_COVARIANCE_TYPES , PN)  , title = ""))
end)
"""

# ╔═╡ 131ceae5-9327-4822-a1d3-f3e2c24ed4c5
md""" 

##### Covariance parameter
$(@bind tau  if !is_individual_settings
	 confirm(NumberField(1e-4:1e-4:10 ,  default = 0.3))
else
	confirm(PF.multi_values_table(PlutoUI.NumberField , projects_names , default_values =rep_tup(1e-4:1e-4:10 , PN) , defaults = rep_tup(0.3 , PN)  , title = ""))
end)
"""

# ╔═╡ c4bb1772-36d2-4963-be39-f53c8e95a29b
md""" 

##### Covariance type
$(@bind reg_type if !is_individual_settings
	 Select(IHT.ALL_REGULARIZATION_TYPES)
else
	confirm(PF.multi_values_table(PlutoUI.Select , projects_names , default_values =rep_tup(IHT.ALL_REGULARIZATION_TYPES , PN)  , title = ""))
end)
"""

# ╔═╡ f3e24d8a-e35d-41de-92e6-0df290d3503c
begin  
	cargs = []
	let (covariance_type , tau) = is_individual_settings ? (covariance_type , tau) : (rep_tup(covariance_type , PN) , rep_tup(tau , PN))
		
		for (ct , tau) in zip(covariance_type , tau)
			if ct <: IHT.AR1Covariance
				push!(cargs , (tau , 1.0))
			elseif ct <: IHT.RelativeDiagonalCovariance
				push!(cargs, (tau, 1.0))
			else
				push!(cargs, Tuple([]))
			end
		end
	end
end;

# ╔═╡ 1318f293-f606-4a9f-8896-853b87c0665d
begin
	xpoints_number = dp_xpoints
	tpoints_number = dp_tpoints
	if is_data_ready
		probls = []
		
		r_types= is_individual_settings ? reg_type : rep_tup(reg_type , PN)
		
		x_p_tuple = is_individual_settings ? xpoints_number : rep_tup(xpoints_number , PN)
		
		t_p_tuple = is_individual_settings ? tpoints_number : rep_tup(tpoints_number , PN)
		
		cov_type_tuple = is_individual_settings ? covariance_type : rep_tup(covariance_type , PN)
		
		for ((_ , d_i), cargs , reg_type , xpoints_number , tpoints_number , covariance_type) in zip(all_data.d, cargs , r_types , x_p_tuple , t_p_tuple , cov_type_tuple)
			
			cov = covariance_type(cargs...)
			inv_probl = IHT.SingleInverseProblem(d_i, C,λ, dλdT, xpoints_number , tpoints_number , covariance=cov , regularization = reg_type())
			push!(probls , inv_probl)
		end
		parallel_probls =IHT.ParallelInverseProblems(Tuple(probls)...)
	end
end;

# ╔═╡ fe228354-5f3f-433d-9f8b-7d888e9c5ecb
begin 
	IHT.set_regularization_multiplier!(parallel_probls , reg_multiplier)
end

# ╔═╡ 6cd83678-860c-41c8-b3cd-007725f9e01d
begin 
	
	bc_plot = Plots.plot(;plot_common_args...)

	for (p_n , inv_probl) in enumerate(parallel_probls.problems)
		t_dir = OHT.trange(inv_probl.direct_problem)
		Plots.plot!(bc_plot, t_dir ,inv_probl.direct_problem.bc_dwn.(t_dir), label = "P$(p_n):lower bc"; plot_common_args...)
		
		Plots.plot!(bc_plot, t_dir ,inv_probl.direct_problem.bc_up.(t_dir), label = "P$(p_n):upper bc"; plot_common_args...)
		for (i,c) in enumerate(eachcol(inv_probl.Tdata_measured))
			Plots.plot!(bc_plot, t_dir, c;plot_common_args..., label="P$(p_n):T$(i)" )
		end
		xlabel!(bc_plot,"Time,s")
		ylabel!(bc_plot,"Temperature, °C")
	end
end

# ╔═╡ 544201c6-90b3-48dd-9651-69ca3c5c4979
bc_plot

# ╔═╡ e9e9e16d-0b2c-45ce-aa5b-cdcda6b143f1
begin 
	reg_multiplier
	is_data_ready
	
	if is_fit_on 
		(start, lb, ub)  =  IHT.fill_starting_vectors(parallel_probls)
		opt_fun = OptimizationFunction(IHT.discrepancy! , 
									  grad = IHT.fdif_gradient!
					)
		optp = OptimizationProblem(opt_fun, start, parallel_probls, lb = lb, ub = ub , maxiters=pso_iters )
		res = solve(optp, optimizer())
		disc_before = IHT.discrepancy!(res.u, parallel_probls)
		if is_after_fit
			new_start = res.u
			disc_before = IHT.discrepancy!(new_start, parallel_probls)
			optp2= OptimizationProblem(opt_fun, new_start, parallel_probls , lb = lb, ub = ub , maxiters=pso_iters)
			res = solve(optp2, NLopt.LN_BOBYQA())
			disc_after = IHT.discrepancy!(res.u, parallel_probls)
		else
			disc_after = 0
			
		end
	else
		res = nothing
	end
	
	
	refresh_graph = false
end

# ╔═╡ 95dbc55f-ab5a-4828-a1e2-9a0c9a9ec19b
begin 
	refresh_graph
	loss_table = IHT.loss_distribution_matrix(parallel_probls)

	md"""
	#### Loss function distribution over problems and types:
	"""
	
end

# ╔═╡ 672119a2-7a47-4813-a2d3-e0c15ee63491
begin 
	pretty_table(HTML,hcat(table_data[:,1] , [v for v in loss_table]...) , column_labels  =vcat("name", [String(k) for k in keys(loss_table)]...) , title = "Loss distribution")
end

# ╔═╡ f0380d5c-ba65-471c-beea-65d766634bc2
if @isdefined res 
	refresh_graph
	result_table(res)
else
	nothing
end

# ╔═╡ 4c9a5a5f-5839-4211-b561-7539dfa74a7a
begin 
	refresh_graph
	
	IHT.loss_distribution(parallel_probls)
end

# ╔═╡ 0c2cf1d2-cb4d-4efa-881c-434784a92366
begin 
	
	refresh_graph
	if is_fit_on && (IHT.optimizable_parnumber(parallel_probls) > 0)
		stats = IHT.ip_covariance(parallel_probls , res.u)
		auto_cor_stats_weighted = IHT.autocorrelation_analysis(parallel_probls)
		auto_cor_stats_unweighted = IHT.autocorrelation_analysis(parallel_probls , is_unweighted = true)
		simple_stats_table = descriptive_stats_table(parallel_probls)
		reg_result_table = regression_stats_table(stats , res.u)

		sensitivity_table = sensitivity_stats_table(parallel_probls , res.u)
	else
		stats = nothing
		auto_cor_stats_weighted = nothing
		simple_stats_table = nothing 
		reg_result_table = nothing
		sensitivity_table = nothing 
	end
	(autocor_weighted_stats_table , autocor_unweighted_stats_table)  = if auto_cor_stats_weighted != nothing
		(
		   autocorrelation_stats_table(auto_cor_stats_weighted),
		 	autocorrelation_stats_table(auto_cor_stats_unweighted , is_unweighted=true)
		)
	else
		(nothing , nothing)
	end

end;

# ╔═╡ 538487b4-b2fc-42c2-ba69-663b2ca5b768
begin 
	c_plot = Plots.plot(Tplot,C.(Tplot) ./density, label = "Inverse", xlabel = "Temperature, °C",linecolor = :green , ylabel = "Heat capacity, J/(kg⋅°C)"; plot_common_args...)
	
	if show_passport
		Plots.plot!(c_plot ,c_x ,c_y ./density , label = "Reference",marker = :circle , markersize  = 8; plot_common_args... )
	end
	if !isnothing(stats) && is_show_cp_conf && is_optimize_c && !is_optimize_lambda

			_o = eval_conf_bounds(stats ,C_poly , α = 2.0)
			plot!(c_plot , Tplot , _o.up.(Tplot)/density , label = nothing)
			plot!(c_plot , Tplot , _o.low.(Tplot)/density , fillrange = (_o.up.(Tplot)/density , _o.low.(Tplot)/density) , alpha = 0.3 , label = "±Δ")
	end
	c_plot
end

# ╔═╡ b81ed1ca-cf5d-40e3-8eaf-903d720ef43e
simple_stats_table

# ╔═╡ eaa0040f-ca36-420b-ae74-d557806114c6
reg_result_table

# ╔═╡ c7cc981a-b8bd-495f-a214-8696de49ee10
autocor_weighted_stats_table

# ╔═╡ ecb7cb75-3257-448a-8bcd-7d4160343fdc
autocor_unweighted_stats_table

# ╔═╡ 8cbe7294-bc96-4386-bf8c-205aa4633f4d
sensitivity_table

# ╔═╡ f2943cf8-ffb9-4065-a272-1f344488dd0f
begin 
	refresh_graph
	problems_number = length(parallel_probls.problems)
	# plotting fitted values 

	p_fit_lam = Plots.plot(Tplot, λ.(Tplot), label = "Inverse", linewidth = 2,linecolor = :green)
	if haskey(L_Dict, material_name) && show_passport
		d_cur = L_Dict[material_name]
		f_range =  (d_cur.x .>= Tmin) .& (d_cur.x .<=Tmax)
		Plots.plot!(p_fit_lam, d_cur.x[f_range],d_cur.y[f_range], label = "Reference",marker = :circle;  plot_common_args...)
	end
	if !isnothing(stats) && is_show_lambda_conf && is_optimize_lambda && !is_optimize_c

			_ol = eval_conf_bounds(stats ,λ.p , α = 1.0)
			plot!(p_fit_lam , Tplot , 
				  _ol.up.(Tplot) , label = nothing)
			plot!(p_fit_lam , Tplot , 
				  _ol.low.(Tplot) , 
				  fillrange = (_ol.up.(Tplot) , 
							   _ol.low.(Tplot)
							  ) , 
				  alpha = 0.3 , 
				  label = "±Δ")
	end
	
	xlabel!(p_fit_lam, "Temperature, ᵒC")
	ylabel!(p_fit_lam, "Thermal conductivity, W/(m*K)")
	#title!(p_fit_lam, "Fitted VS measured")

	# residuals
	
	p_distr = Plots.plot(;plot_common_args...)

	is_show_weighted_residuals = (residuals_type == :weighted) && (residuals_type != :both)
	is_show_both = residuals_type == :both

	_N = !is_show_both ? 1 : 2
	
	p_residual_tuple = ntuple(_N) do _
			Plots.plot(;plot_common_args...)
	end
	p_autocor_tuple = ntuple(_N) do _
		Plots.plot(;plot_common_args...)
	end
	p_hist_tuple = ntuple(_N) do _
		Plots.histogram(;plot_common_args...)
	end
	p_jac_plot_tuple = ntuple(length(parallel_probls.problems)) do _
			Plots.plot(;plot_common_args... , xlabel = "timestep")
	end
	_title = if is_show_weighted_residuals
		(("Weighted residuals" , 
		 	"Weighted residuals",
		 	"Weighted residuals autocorrelation"),)
	elseif is_show_both
		(
			("Weighted residuals" , "Weighted residuals" , "Weighted residuals autocorrelation") ,
			("Residuals" , "Residuals" , "Residuals autocorrelation")  
		 	  
		)
	else
		(
			("Residuals" , 
		 	 	"Residuals" , 
		 		"Residuals autocorrelation"),
		)
	end
	
	for (p_n , inv_probl) in enumerate(parallel_probls.problems)
		
		# plotting temperature disctribution
		for (i,(te,tm)) in enumerate(zip(eachcol(inv_probl.Tdata_evaluated) , eachcol(inv_probl.Tdata_measured)))
			Plots.plot!(p_distr , te, label = "P$(p_n):T$(i)clc")
			Plots.plot!(p_distr, tm, linestyle = :dash, label = "P$(p_n):T$(i)exp"; plot_common_args...)
		end

		title!(p_distr , "Temperature distributions")
		xlabel!(p_distr, "Timestep index")
		ylabel!(p_distr, "Temperature, ᵒC")
		# plotting jacobians 
		local_minima = !isnothing(res) ?  res.u : IHT.fill_starting_vectors(parallel_probls).x₀
		if length(local_minima) > 0 
			J = IHT.fdif_jacobian(inv_probl , local_minima) 
			for (i , c) in enumerate(eachcol(J))
				plot!(p_jac_plot_tuple[p_n] , c,  label = "β$(i)" )
			end
		end
		title!(p_jac_plot_tuple[p_n] , "Problem $(p_n)" )
			
		# plotting residuals
		res_iterator = if is_show_weighted_residuals
			res_iterator = enumerate(zip(
				eachcol(IHT.ResidualIterator(Val(true), inv_probl)), 
			)
			)
		elseif is_show_both
			enumerate(
				zip(
					eachcol(IHT.ResidualIterator(Val(true), inv_probl)) , eachcol(inv_probl.residual)
				)
			)
		else
			enumerate(zip(
				eachcol(inv_probl.residual),
								)
					 )
		end

		for (i , c) in res_iterator
			
			
			cur_label = "P$(p_n):T$(i)"
			pargs = (label = cur_label,)
			
			for (c_ii, p_r , p_a , p_h , ttl) in zip(c , 
											   p_residual_tuple ,
											   p_autocor_tuple , 
											   p_hist_tuple , 
												 _title)
				
				c_i = collect(c_ii)
				Plots.plot!(p_r ,  c_i; pargs... , title = ttl[1])
		
				lgs = collect(0:(length(c_i) - 1))
				
				Plots.plot!(p_a , autocor(c_i , lgs) ; pargs..., title = ttl[3])
	
				bb = range(minimum(c_i),maximum(c_i),20)
				Plots.histogram!(p_h , c_i ; bins = bb,						 alpha=0.5 , pargs..., title = ttl[2])
			end
			
		end
	end
	if is_show_both
		_args = (plot_common_args... , layout = (2,1) )
		p_residual = Plots.plot(p_residual_tuple...; _args...)
		p_autocor =  Plots.plot(p_autocor_tuple...; _args...)
		p_hist =  Plots.plot(p_hist_tuple...; _args...)
	else
		p_residual = first(p_residual_tuple)
		p_autocor = first(p_autocor_tuple)
		p_hist = first(p_hist_tuple)
	end
	Plots.plot!(p_distr;plot_common_args...)
	  Plots.plot!(p_autocor;plot_common_args...)
	  Plots.histogram!(p_hist;plot_common_args...)

	xlabel!(p_residual, "Timestep")
	ylabel!(p_residual, "ΔT, ᵒC")
	ylabel!(p_hist , "Counts number")
	xlabel!(p_hist , "ΔT")
	xlabel!(p_autocor, "Timestep")

	p_jac_plot = Plots.plot(p_jac_plot_tuple... ;plot_common_args... , layout = (problems_number , 1))

end ;

# ╔═╡ 67763d8f-e1ac-4b1f-978f-4734af1e03ba
ylims!(p_fit_lam, lam_y_scale_region)

# ╔═╡ 33280e2b-2567-4289-970a-b034ba4c8cd2
md""" 
### Residuals Histogram
--------------
$p_hist
"""

# ╔═╡ a660bf26-5b50-4910-b0ab-8e453623dc1a
md""" 
### Residuals time trend
--------------
$p_residual
"""

# ╔═╡ 9464897c-8591-4fe3-aa24-9c710a37f7b6
md""" 
### Residuals autocorrelation function
--------------
$p_autocor
"""

# ╔═╡ 042ea29b-63dd-43d0-a20f-68807c5f7cd4
md""" 
### Temperature trends
--------------
$p_distr
"""

# ╔═╡ 5f6f0710-4943-4fb4-b56a-251f74ce584d
md""" 
### Sensitivity plots
--------------
$p_jac_plot
"""

# ╔═╡ dd8fb4fa-fa67-4e26-988d-3ede79bc9540
if use_l_curve
	#res_l = is_after_fit ? res2 : res 
	
	alpha_range = collect(range(extrema(l_curve_alpha)... , length=l_curve_points_number))
	
	alphas = is_l_curve_logscale ? 10.0 .^alpha_range : alpha_range

	lcurve_probs = IHT.ParallelInverseProblems(Tuple(probls)...)
	
	(_start, _lb, _ub)  =  IHT.fill_starting_vectors(lcurve_probs)
	
	
	sols = Dict{Int , Any}()
	Threads.@sync for i in 1:length(alphas)
		 Threads.@spawn begin
			p_i = deepcopy(lcurve_probs)
			IHT.set_regularization_multiplier!(p_i , alphas[i])
			optp_i = OptimizationProblem(IHT.discrepancy!, 
										 _start, 
										 p_i,
										 lb = _lb, 
										 ub = _ub)
			sols[i] = solve(optp_i, optimizer()  , maxiters=pso_iters )
		end
	end

	
end

# ╔═╡ cf8665e2-07d9-4464-b833-0531205408e9
if  use_l_curve
	l_curve_cov_loss = similar(alphas)
	l_curve_reg_loss = similar(alphas)
	l_curve_total_loss = similar(alphas)
	for (i , r) in sols
		IHT.discrepancy!(r.u , lcurve_probs)
		_loss_full = IHT.loss_distribution(lcurve_probs)
		l_curve_cov_loss[i] = _loss_full.covariance
		l_curve_reg_loss[i] = _loss_full.regularization/alphas[i]
		l_curve_total_loss[i] = _loss_full.total
	end
	l_curve_ind = sortperm(l_curve_cov_loss)
	permute!(l_curve_cov_loss, l_curve_ind)
	permute!(l_curve_reg_loss,l_curve_ind)
	permute!(l_curve_total_loss, l_curve_ind)
	permute!(alphas,l_curve_ind)
end

# ╔═╡ 2614e5ad-0d45-4046-abfd-053dd872c17c
if use_l_curve
		best_alf_Index = Main.find_l_corner_index(log10.(l_curve_cov_loss) , log10.(l_curve_reg_loss))

	l_curve_best_alpha = alphas[best_alf_Index]

end

# ╔═╡ 3f8c5c6a-f902-4b1e-abc1-7b572c8d2505
if use_l_curve
	md"""
	##### L-curve methods suggests to try α = $(l_curve_best_alpha)
	"""
end

# ╔═╡ 1c32440c-32af-414c-9a7a-d6cc83685f7c
if  use_l_curve 
	lessem16 = t -> t <= 1e16
	scale_type = (any(lessem16 , l_curve_cov_loss) || any(lessem16 , l_curve_reg_loss) ) ? :identity : :log10

	scale_args = (xscale = scale_type, yscale = scale_type)
	

	
	reg_plot = Plots.scatter(;plot_common_args...)
	Plots.scatter!(reg_plot , alphas , l_curve_cov_loss , label = "model " ; scale_args... , plot_common_args...)
	
	Plots.scatter!(reg_plot , alphas , l_curve_reg_loss , label = "reg " ; scale_args..., plot_common_args...)
	
	Plots.scatter!(reg_plot , alphas , l_curve_total_loss, label = "total " ; scale_args..., plot_common_args...)

	
	min_min_y = min(minimum(l_curve_best_alpha), minimum(l_curve_cov_loss))
	max_max_y = maximum(l_curve_total_loss)

	Plots.plot!(reg_plot ,[l_curve_best_alpha , l_curve_best_alpha]  , [min_min_y , max_max_y] , label = "best" , legend_position = :best ; scale_args..., plot_common_args...)

	
	xlabel!(reg_plot , "Regularization multiplier , α")
	ylabel!(reg_plot , "log(loss)")
	title!(reg_plot , " Best α = $(l_curve_best_alpha)")
	
	l_curve_plot = Plots.plot(;plot_common_args...)
	Plots.scatter!(l_curve_plot , l_curve_cov_loss , l_curve_reg_loss ; scale_args..., plot_common_args... , label = nothing)
	xlabel!(l_curve_plot , "Model prediction loss")
	ylabel!(l_curve_plot , "Regularization loss"  )
	
end;

# ╔═╡ 90c287a9-5719-4495-9a3a-c5ac21817cb3
use_l_curve && l_curve_plot

# ╔═╡ ccbd5672-7385-4740-8e9f-d46a7afa8753
use_l_curve && reg_plot

# ╔═╡ c2564f33-dcf5-45fe-a462-70fe105bcf0a
if false 
	prob_fun =  IHT.regularization_scan(alphas , Optimization.remake)
	opt_l_curve_probl = OptimizationProblem(IHT.discrepancy!, 
										 _start, 
										 deepcopy(lcurve_probs),
										 lb = _lb, ub = _ub)
	ensmbl =  Optimization.EnsembleProblem(
                    opt_l_curve_probl, 
                    prob_func = prob_fun
    ) 

	res_ensemble = solve(ensmbl, optimizer(), EnsembleThreads(), trajectories = length(alphas) ,  maxiters=pso_iters)

end

# ╔═╡ f5b13ec5-98f9-48bb-ad95-7dbd87e11f7b
lam_test = OHT.thermal_conductivity(parallel_probls.problems[1].direct_problem, Ttest)

# ╔═╡ 27031733-7ade-4bda-b464-5a9ade0b2950
begin 
	name_poly =  material_name
	x_fit = L_Dict[name_poly].x
	www =   (d_cur.x .>= Tmin) .& (d_cur.x .<=Tmax)
	x_fit = x_fit[www]
	y_fit = L_Dict[name_poly].y[www]
	(k,_,rrr) =  PW.scale_x_to_ξ(x_fit)
	pp = PW.polyfit(PW.BernsteinSymPoly{lam_basis_degree, Float64} , k , y_fit)
	ppp = Plots.plot(x_fit,pp.(k), label = nothing)
	Plots.plot!(ppp,x_fit,y_fit, marker=:circle, label = nothing; plot_common_args...)
	Plots.scatter!([Ttest] ,[lam_test] )
	title!(ppp,"""Polynomial fitting of passport,
		   """)
end

# ╔═╡ 4cf77d64-a720-41da-a9c2-5d875ea00135
if is_write_file 
	resave_file
	refresh_graph
	_ff_init = joinpath(data_selection_path_ref[] , data_selection_name_ref[])
	ff_name = joinpath(save_folder , default_save_file)
	cp(_ff_init, ff_name; force = true)
	DC.WinPos.export_to_hdf5(parallel_probls, ff_name)
	"✅ Data selection saved to hdf5-file $(ff_name) at $(Dates.format(now(), "HH:MM:SS"))"
end

# ╔═╡ ea619c49-c928-494d-827f-185b9137f27e
IHT.coeffs(parallel_probls.problems[1].optimizable.λ)

# ╔═╡ 2c692ea0-c6d7-4872-9d8d-f1c501e734ec
parallel_probls

# ╔═╡ Cell order:
# ╟─a17fe1fe-5542-454b-b45e-942ac52b6f1a
# ╟─5807712b-5d26-49c8-ab65-dac167ebad7b
# ╟─2bfb4e52-6248-4832-aca1-98ba58959bff
# ╟─3b1c3b0a-558e-4987-bf16-072963e455cf
# ╟─2db7cd7e-3d1a-4533-b0de-27188300565b
# ╟─db671921-13dc-497b-81e5-dcb4da0695f9
# ╟─450fb200-eec6-4e96-9ebd-81453c015830
# ╟─6e062bd9-d20c-4e1d-b772-328bec8859ea
# ╟─a7abe643-2553-450d-ac81-d4a690a1c2ff
# ╟─41bc1a0a-73c8-430d-a1d3-4eb98487c815
# ╟─34b63c47-678d-4a5f-aed3-1896b778e117
# ╟─444236e5-e010-4b8b-8709-31c0307dc5d8
# ╟─81c2f93b-05b1-4eb0-9919-4ef76ecad233
# ╟─074c47c6-a1ae-4ded-8e64-b40393d7ba4a
# ╟─35b5c27e-921f-4fe7-90bc-3b327d9150fe
# ╟─5be15388-cea2-4884-b8de-bff5be64e506
# ╟─6c4cb363-3bda-4dfa-8499-8201a013895d
# ╟─b31a7c52-7c4b-480f-84d8-fcb1c71dca1b
# ╟─c69d07a3-ba0a-48a2-b296-1944b8cd322c
# ╟─b1f00c55-5ce9-4a0f-a548-a8a9041d02fc
# ╟─544201c6-90b3-48dd-9651-69ca3c5c4979
# ╟─62069546-44fb-4a77-986d-f03624719e29
# ╟─8f65e989-1abe-4689-9c3f-fdb6cabd7eae
# ╟─a46970c8-7c91-4795-83a9-41c56a7ca399
# ╟─7d37019a-34bf-43ca-aa81-8b6c29191473
# ╟─74707b89-d91c-468d-9b4a-a06dc99f69c0
# ╟─4f3e1899-541d-4809-810c-f2e6b7ca1ad1
# ╟─4fdbfb71-ec8d-4d30-b7e7-f289f773fadd
# ╟─eb755603-1eda-4182-a026-e9b163e78ae3
# ╟─97618f0c-52ae-4772-b5b7-5512ea44af09
# ╟─baffb68a-fb52-4bac-b484-4a215108aaed
# ╟─c8cda804-9a2a-40c9-aa91-7a48500e85ff
# ╟─538487b4-b2fc-42c2-ba69-663b2ca5b768
# ╟─8345302e-0b9d-4208-9a66-3d9f32903b39
# ╟─5843fcfb-4e0e-480d-b263-e3f3ad6a7ac3
# ╟─d051f5e5-f836-4043-9326-639b40acaf87
# ╟─16337bb4-438e-4b86-bdc5-b88ef210a960
# ╟─0ff89750-2a97-436b-b759-7da1354b2c6f
# ╟─4ca95124-8a7a-4e4f-9e65-ff7b2adf35a5
# ╟─fa72774e-040a-4bc3-a759-5eb68c243fb4
# ╟─c8dc4f9d-a549-4dc2-82bd-38ffe949ea55
# ╟─bfa23359-8bac-4db0-bac1-0885ebe8ec4b
# ╟─88e8d37d-e4e4-486d-921e-03a74fbf00f2
# ╟─67763d8f-e1ac-4b1f-978f-4734af1e03ba
# ╟─68ed85b2-7308-4ea7-b696-0f1951219592
# ╟─4330cdcd-24fc-459b-8d29-a935c6a7c347
# ╟─ce3dc022-0f15-41cb-8cd2-2c33b726c482
# ╟─c8861fa2-9cb3-4349-9ece-eba285c24eab
# ╟─9b35c13b-24ba-4c43-a621-a3f8ba45fe4a
# ╟─fc5c1209-26ec-41d7-9238-e57d18330de1
# ╟─2b3a1a65-8d3c-424e-a6cf-0e96646795f4
# ╟─3432fabd-911b-4370-af86-c396cbc7bbab
# ╟─3281dd3c-0453-4098-a6ed-4d771717033b
# ╟─63f66c3f-1489-4da1-963b-47812ca245ae
# ╟─95dbc55f-ab5a-4828-a1e2-9a0c9a9ec19b
# ╟─107318d8-051f-4967-8dc0-49d836cd9469
# ╟─131ceae5-9327-4822-a1d3-f3e2c24ed4c5
# ╟─c4bb1772-36d2-4963-be39-f53c8e95a29b
# ╟─33473b7a-e22f-4333-a2f2-374778c0d603
# ╟─aad5d954-fb91-4a7d-a09e-02105ef3d0f9
# ╟─672119a2-7a47-4813-a2d3-e0c15ee63491
# ╟─b81ed1ca-cf5d-40e3-8eaf-903d720ef43e
# ╟─eaa0040f-ca36-420b-ae74-d557806114c6
# ╟─c7cc981a-b8bd-495f-a214-8696de49ee10
# ╟─ecb7cb75-3257-448a-8bcd-7d4160343fdc
# ╟─5e5b27d1-fdc8-451d-a4f0-93f33adbcf76
# ╟─33280e2b-2567-4289-970a-b034ba4c8cd2
# ╟─a660bf26-5b50-4910-b0ab-8e453623dc1a
# ╟─9464897c-8591-4fe3-aa24-9c710a37f7b6
# ╟─042ea29b-63dd-43d0-a20f-68807c5f7cd4
# ╟─5f6f0710-4943-4fb4-b56a-251f74ce584d
# ╟─8cbe7294-bc96-4386-bf8c-205aa4633f4d
# ╟─5d4b49c0-1d0f-41b5-b359-30c033fc8544
# ╟─2a1ca349-3732-443e-bc48-5be611a5d91f
# ╟─d02ec3fd-1b0d-4bc3-92e9-adedc8bf2c8c
# ╟─0a0324df-430f-4ac0-857b-4da6a7dca138
# ╟─f0380d5c-ba65-471c-beea-65d766634bc2
# ╟─c8cd1797-b631-4fdf-a51f-67e22c86c55a
# ╟─aea6fc8f-5ee2-40a3-aebd-942c45eec0d6
# ╟─eb70f6d9-301f-4e0f-91a2-0e911b1a6dd9
# ╟─a8f5cc07-b2bd-445c-8797-ffe78a25641b
# ╟─10f7c083-a697-4b7e-87ea-de2791ed1a30
# ╟─3f8c5c6a-f902-4b1e-abc1-7b572c8d2505
# ╟─90c287a9-5719-4495-9a3a-c5ac21817cb3
# ╟─ccbd5672-7385-4740-8e9f-d46a7afa8753
# ╟─fe228354-5f3f-433d-9f8b-7d888e9c5ecb
# ╟─f3e24d8a-e35d-41de-92e6-0df290d3503c
# ╟─4c9a5a5f-5839-4211-b561-7539dfa74a7a
# ╟─fc79d398-03d0-4f75-a777-41e2e27fee5e
# ╠═1318f293-f606-4a9f-8896-853b87c0665d
# ╟─6cd83678-860c-41c8-b3cd-007725f9e01d
# ╟─f2943cf8-ffb9-4065-a272-1f344488dd0f
# ╠═ffeafacd-2b98-48c0-b840-297fb51f5e54
# ╟─ea412a80-0bf7-4ac6-9b90-8530a4c26008
# ╟─23902c7d-5973-4868-8488-e0c7634573c4
# ╟─6800ae42-c5b1-4d1e-84fb-a03015bf138f
# ╟─c4053281-7cd4-4590-b004-b4ca43771094
# ╠═9cd8c4c8-7d11-4fb0-92d5-439702aa9496
# ╠═e9e9e16d-0b2c-45ce-aa5b-cdcda6b143f1
# ╠═0c2cf1d2-cb4d-4efa-881c-434784a92366
# ╟─70e2c8f2-d896-4773-8280-d391d9975307
# ╟─bce5f7ae-49c5-4520-a0f6-6215e5078674
# ╟─dd8fb4fa-fa67-4e26-988d-3ede79bc9540
# ╟─cf8665e2-07d9-4464-b833-0531205408e9
# ╟─2614e5ad-0d45-4046-abfd-053dd872c17c
# ╟─1c32440c-32af-414c-9a7a-d6cc83685f7c
# ╟─c2564f33-dcf5-45fe-a462-70fe105bcf0a
# ╟─27031733-7ade-4bda-b464-5a9ade0b2950
# ╟─f5b13ec5-98f9-48bb-ad95-7dbd87e11f7b
# ╟─d71fd6d1-167d-40fb-a253-b0982e19c0d0
# ╟─4cf77d64-a720-41da-a9c2-5d875ea00135
# ╟─58e0b175-a114-46be-9792-35954faea5ec
# ╟─efa9120f-5a45-41a0-9132-dc26f967fec3
# ╟─fb2937d0-738b-4329-aef5-f3af1d17497a
# ╟─023dc25f-6cf8-4802-83b4-77d1827cd2a7
# ╟─ab211f5e-9cf8-438c-a3cc-4fa18960b67b
# ╟─2bf91a7e-0da8-48e8-a779-962be2e7c03b
# ╟─fd51a6c8-6569-4bfd-86ac-883c648fe6d9
# ╟─dbd6f9d9-4b44-4238-86ca-7cd361545a56
# ╟─206e1e8f-16f3-4144-8401-c2cbf40c0125
# ╟─4b8ddfc0-ed4f-4b66-b752-fa075d348608
# ╟─04ad22b9-e474-41ed-bc87-7dbf9a2b1dcc
# ╟─3b9e739c-9519-4efa-b2da-cac5451d55d3
# ╟─f799c1f9-ebc2-4fdb-839e-ee721ab9b8f5
# ╠═544c4618-5b5e-4b25-a916-a43fc2e05913
# ╠═e83d394a-f06f-48de-a664-aad82e4361fc
# ╠═baaf2f95-9529-4829-a8fa-0724fd21b78e
# ╠═b32b35a1-f31c-428f-a9ab-d640f2be496e
# ╟─9c4936fd-25d6-4760-85f5-38aba54f1800
# ╠═ea619c49-c928-494d-827f-185b9137f27e
# ╠═2c692ea0-c6d7-4872-9d8d-f1c501e734ec
