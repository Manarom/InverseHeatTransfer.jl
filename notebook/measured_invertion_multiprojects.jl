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

# ╔═╡ a17fe1fe-5542-454b-b45e-942ac52b6f1a
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
end

# ╔═╡ 2bb45200-d7ba-47a3-b82d-9c147b5d7601
default_data_fodler = joinpath(@__DIR__, "..","test","test_data","property_inversion_ansys_new")

# ╔═╡ 80479ac3-897a-4aa7-8ce9-977fa4912bc6
source_path = joinpath(@__DIR__,"..","src")

# ╔═╡ 6d2a48d4-e458-4edd-b22a-7b7dae9f492c
#includet(joinpath(source_path, "WinPos.jl"))

# ╔═╡ a6d591eb-20f2-476f-af19-9b0084ddf929
IHT = InverseHeatTransfer

# ╔═╡ 35958e8a-eb7a-4eff-89a0-f9c04aff2a37
begin 
	
	PW = IHT.ScaledPolynomials
	OHT = IHT.OneDHeatTransfer
	WP = IHT.WinPos
	DC = IHT.DataConnector
end

# ╔═╡ 438a3909-9367-4660-a3ca-bd1786ab6016
plot_common_args = (grid = true, gridlinewidth=3, gridstyle = :dot,minorgrid=true, box = :on, linewidth = 3)

# ╔═╡ db671921-13dc-497b-81e5-dcb4da0695f9
md""" ## Loading experimental data from winpos project"""

# ╔═╡ 450fb200-eec6-4e96-9ebd-81453c015830
md" Load data from : $(@bind input_data_type Select([:winpos, :ascii , :hdf5] , default = :files))"

# ╔═╡ a407a99b-b40c-436c-a2a0-af2e19b347b8
is_winpos = input_data_type == :winpos;

# ╔═╡ 6e062bd9-d20c-4e1d-b772-328bec8859ea
begin 
	md""" #### working folder $(@bind working_folder TextField(90, default = realpath(default_data_fodler))) """
end

# ╔═╡ 81c2f93b-05b1-4eb0-9919-4ef76ecad233
begin 
	projects = WP.load_from_winpos_folder(working_folder)
	#projs_txt  = WP.load_from_ascii_folder(working_folder ; name_matcher="Tmeasured", variable_name = "T")
	#projects = merge(projs_wp , projs_txt)
end;

# ╔═╡ 17fbc55f-12a8-431e-ac00-b28304f2eb6c
md""" #### select projects to be taken into account: 


$(@bind cur_proj confirm(MultiSelect(collect(keys(projects)) , default = ["10ks" , "2ks"]))) 


"""

# ╔═╡ 0a4ce046-5b43-425d-a3f2-da8d9867b181
#=mutable struct DataConnector
	all_names::Vector{String}
	selected_names::Vector{String}
	data :: Matrix{Float64}
	data_cutted :: Matrix{Float64}
	locations :: Vector{Float64}
	total_thickness :: Float64
	project:: WP.WinPosProject
	
end=#

# ╔═╡ dcf0034b-e405-4f27-b853-acb7a58b9cc6
t_cutted_min(d::DataConnector) = minimum(d.data_cutted[:,1])

# ╔═╡ 4d27d75e-5e95-43c9-9e77-b9ae3d6c4bb9
t_cutted_max(d::DataConnector) = maximum(d.data_cutted[:,1])

# ╔═╡ 28b91e49-996c-4abc-9db4-30b515444aab
function cut_data!(d::DataConnector , tmin , tmax)
	t = @view d.data[: , 1]
	f =@. (t < tmax) & (t > tmin) 
	d.data_cutted = d.data[f,:]
end

# ╔═╡ 054d932d-12be-4538-ab34-b5d3f465bf0f
projects

# ╔═╡ 2e676f8c-909a-4fb7-b35e-57d08f802df7
data_changed = true

# ╔═╡ 6c4cb363-3bda-4dfa-8499-8201a013895d
@bind show_data_table Select(collect(keys(all_data)))

# ╔═╡ fd47bc58-8a0c-4dd7-875c-bcc80a21e64e
all_data

# ╔═╡ 62069546-44fb-4a77-986d-f03624719e29
md" ## Physical properties setup"

# ╔═╡ fc5c1209-26ec-41d7-9238-e57d18330de1
@bind refit Button("Refit!!")

# ╔═╡ 2b3a1a65-8d3c-424e-a6cf-0e96646795f4
md" ### refit $(@bind is_fit_on CheckBox(default = false))"

# ╔═╡ a46970c8-7c91-4795-83a9-41c56a7ca399
md""" ### Compare to passport $(@bind show_passport CheckBox(default = true))"""

# ╔═╡ 4f3e1899-541d-4809-810c-f2e6b7ca1ad1
md"#### Optimize C : $(@bind is_optimize_c CheckBox(default = false))"

# ╔═╡ 2bb61e43-384c-48ce-8a55-843726bf3f05
md" Upper Cp limit = $(@bind upper_cp_limit confirm(Slider(100.0 : 1e-1 : 10000.0, default = 1500.0, show_value = true)))"

# ╔═╡ baffb68a-fb52-4bac-b484-4a215108aaed
md" Lower Cp limit = $(@bind lower_cp_limit confirm(Slider(100.0 : 1e-1 : 10000.0, default = 600.0, show_value = true)))"

# ╔═╡ 4ca95124-8a7a-4e4f-9e65-ff7b2adf35a5
md" #### Optimize λ ? $( @bind is_optimize_lambda CheckBox(default = true))"

# ╔═╡ c8dc4f9d-a549-4dc2-82bd-38ffe949ea55
md" Lower λ limit = $(@bind lower_lam_limit confirm(Slider(0.0 : 1e-3 : 10.0, default = 0.1, show_value = true)))"

# ╔═╡ bfa23359-8bac-4db0-bac1-0885ebe8ec4b
md" Upper λ limit = $(@bind upper_lam_limit confirm(Slider(0.1 : 1e-3 : 30.0, default = 20, show_value = true)))"

# ╔═╡ 68ed85b2-7308-4ea7-b696-0f1951219592
@bind lam_y_scale_region PlutoUI.combine() do Child 
md"""
	``\lambda_{min} ``= $(
		Child(NumberField(0:1e-2:50,default=0.0))
	) \
	``\lambda_{max} `` = $(
		Child(NumberField(0:1e-3:50,default=30.0))
	)\
	"""
end

# ╔═╡ 9b35c13b-24ba-4c43-a621-a3f8ba45fe4a
md" ### Inverse problem setup"

# ╔═╡ 33473b7a-e22f-4333-a2f2-374778c0d603
md"""

Covariance type $(@bind covariance_type Select(IHT.ALL_COVARIANCE_TYPES))

Covariance parameter τ = $(@bind tau confirm(Slider(1e-3:1e-2:10 , show_value = true, default = 1.0)))

Regularization type : $(@bind reg_type Select(IHT.ALL_REGULARIZATION_TYPES))

Regularization multiplier α: $(@bind reg_multiplier confirm(NumberField(0.0 : 1e-3 : 1000 , default = 1e-3)))
"""

# ╔═╡ c8cd1797-b631-4fdf-a51f-67e22c86c55a
md"""
	### L-curve analysis $(@bind use_l_curve CheckBox(default = false))

	α points number $(@bind l_curve_points_number Select(3:100 , default = 20) )
	"""

# ╔═╡ 0a0324df-430f-4ac0-857b-4da6a7dca138
md"addtitional fit $(@bind is_after_fit CheckBox(default = false))"

# ╔═╡ f3e24d8a-e35d-41de-92e6-0df290d3503c
begin 
	if covariance_type <: IHT.AR1Covariance
		cargs = (tau , 1.0)
	elseif covariance_type <: IHT.RelativeDiagonalCovariance
		cargs  = (tau, 1.0)
	else
		cargs = Tuple([])
	end
end;

# ╔═╡ fa72774e-040a-4bc3-a759-5eb68c243fb4
md" Λ(T) parameters number $(@bind basis_degree Select(1:100, default = 4))"

# ╔═╡ 3281dd3c-0453-4098-a6ed-4d771717033b
md" Direct problem number of coordinate steps $(@bind dp_xpoints Select(10:2000, default = 140))"

# ╔═╡ bbc2d0df-28bf-4754-bbdf-bece0bbfe76b
md" Direct problem number of time steps $(@bind dp_tpoints Select(10:10:5000, default = 2000))"

# ╔═╡ 2a1ca349-3732-443e-bc48-5be611a5d91f
md"""
### Select optimizer $(@bind optimizer Select([ParticleSwarm => "Particle swarm", OptimizationNLopt.NLopt.LN_COBYLA => "COBYLA", OptimizationNLopt.NLopt.LN_BOBYQA => "BOBYQA"  ], default = OptimizationNLopt.NLopt.LN_COBYLA)) 
"""

# ╔═╡ d02ec3fd-1b0d-4bc3-92e9-adedc8bf2c8c
md" ##### Global optimizer iterations number $(@bind pso_iters Select(10:10:10000 , default = 200))"

# ╔═╡ 9cd8c4c8-7d11-4fb0-92d5-439702aa9496
md" ### OPTIMIZATION "

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

#### log10(α) search range: $(@bind l_curve_alpha confirm(RangeSlider(l_culve_alphas_default_range)) )


"""

# ╔═╡ 10f7c083-a697-4b7e-87ea-de2791ed1a30
md" α = $(10 .^extrema(l_curve_alpha)) " 

# ╔═╡ bce5f7ae-49c5-4520-a0f6-6215e5078674
includet(joinpath(source_path, "problem_ensemble_functions.jl"))

# ╔═╡ 84e23b57-e5a4-4bf0-97ef-6b41b502b4e6
optimizer

# ╔═╡ fb2937d0-738b-4329-aef5-f3af1d17497a
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

# ╔═╡ 8f65e989-1abe-4689-9c3f-fdb6cabd7eae
md" #### Select material name : $(@bind material_name Select(collect(keys(Cp_Dict))))"

# ╔═╡ efa9120f-5a45-41a0-9132-dc26f967fec3
Rho_Dict = Dict(
	
	"RBSN" => 2700.0, 
	"NIASIT" => 2000.0,
	"OTM357" => 2530.0,
	"HAFS" => 1600.0,
	"BAFS" => 1800.0
	
);

# ╔═╡ 7d37019a-34bf-43ca-aa81-8b6c29191473
md"""

### Density, kg/m³ $(@bind density confirm(NumberField(200.0:0.1:4000, default = Rho_Dict[material_name])))

"""

# ╔═╡ 5843fcfb-4e0e-480d-b263-e3f3ad6a7ac3
if density != Rho_Dict[material_name]
	md" passport density is $(Rho_Dict[material_name])"
end

# ╔═╡ 023dc25f-6cf8-4802-83b4-77d1827cd2a7
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
)

# ╔═╡ 23902c7d-5973-4868-8488-e0c7634573c4
@bind default_save_file TextField(60,"thermal_conductivity.txt")

# ╔═╡ 6800ae42-c5b1-4d1e-84fb-a03015bf138f
@bind save_folder TextField(90, default = realpath(default_data_fodler))

# ╔═╡ ea412a80-0bf7-4ac6-9b90-8530a4c26008
md"### Save data to file ? $(@bind is_write_file Button())"

# ╔═╡ 2bf91a7e-0da8-48e8-a779-962be2e7c03b
function plot_optimizable(C; xmin,xmax , kwargs...)
	x = range(xmin , xmax ,100)
	return (Plots.plot(x, C.(x); kwargs...), x)
end

# ╔═╡ fba5bc8b-25cb-406b-aa13-f06c591e08c9
function discrepancy(x,inv_probl::IHT.SingleInverseProblem)
	λ = inv_probl.direct_problem.L_f.fun
	IHT.refill!(λ, x)
	dλdT = inv_probl.direct_problem.Ld_f.fun
	IHT.derivative!(dλdT, λ)
	IHT.solve_direct_problem!(inv_probl)
	IHT.fill_residual!(inv_probl)
	#return sum(Base.Fix2(^,2.0), inv_probl.residual)
	return norm(inv_probl.residual)/length(inv_probl.residual) #+ sum(Base.Fix1(^,2.0),dλdT.p.poly.coeffs)/PW.parnumber(dλdT.p)
end

# ╔═╡ 5135015e-0fcf-484a-86ad-3cb7e40849bd
function select_folder_gtk()
    # open_dialog_native или open_dialog с флагом выбора папки
    # В Gtk.jl для папок используется специфический метод
    return open_dialog("Выберите папку", action=GtkFileChooserAction.SELECT_FOLDER)
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

# ╔═╡ 6cd4f554-7242-4e97-b4cb-8549e3b70139
function multi_values(names, default_values=nothing)
	isnothing(default_values) && (default_values = zeros(length(names)))
	PlutoUI.combine() do Child
		@htl("""
		<h6>Thicknesses, mm</h6>
		<ul>
		$([
			@htl("<li>$(name): $(Child(name, NumberField(0:1e3:20, default=deflt)))</li>")
			for (name,deflt) in zip( names, default_values)
		])
		</ul>
		""")
	end
end

# ╔═╡ 646aebc7-b3db-443f-888f-800d444b4fa3
function fill_data_connector_data!(d::DataConnector)
	isempty(d.selected_names) && return false
	(ttt, TTT, winpos_data_names) = Main.WinPos.joindata(d.project, names = d.selected_names)
	d.selected_names = winpos_data_names
	d.data =  hcat(ttt,TTT)
	d.data_cutted = copy(d.data)
end

# ╔═╡ 60458134-1de0-475e-ba64-24b6f33a2980
function paired_selected_names(all_data , field_name :: Symbol = :selected_names)
	out = Vector{Pair{String , String}}()
	for (n,v) in all_data
		append!(out , [ Pair(n , k)  for k in  getfield(v , field_name) ] )

	end
	return out
end

# ╔═╡ 6d0d7cb4-45fc-405f-8a5d-cf10ad5e380b
function multi_values(all_data::AbstractDict, field_name::Symbol =  :selected_names , default_values=nothing)
	PlutoUI.combine() do Child
		@htl("""
		<h6>field_name, mm</h6>
		<ul>
		$([
			@htl("<li>$(join( [n , f] , ":")): $(Child( join( [n , f] , ":"), NumberField(0:1e-3:20 , default = 0.0)))</li>")
			
			for (n, f) in paired_selected_names(all_data)
		])
		</ul>
		""")
	end
end

# ╔═╡ 96b4c7c4-30a8-407c-9bd4-245f0ee0d9b2
function project_selector(all_data , field_name::Symbol)

		PlutoUI.combine() do Child
		@htl("""
		<ul>
		$([
			@htl("<li>$(name): \n $(Child(
				name , 
				MultiSelect(collect(getfield(v , field_name) ))
			)
								   )</li>")
			for (name, v) in all_data
				 ])
		</ul>
		""")
	end

end

# ╔═╡ c30f95e5-93eb-4ab7-bc84-4bfe7092cf47
@bind selected_variables_multi  confirm(project_selector(all_data , :all_names))

# ╔═╡ 0177413c-0f89-4935-9963-5aeebd333b9a
is_selected = !any(isempty , selected_variables_multi)

# ╔═╡ bc6ecff4-2ade-4436-9630-be573eb1ea04
begin 
	is_selected
	@bind thicknesses_mm confirm(multi_values(collect(keys(all_data))))
end

# ╔═╡ ca69ae4b-8428-4643-93a4-c4c1687b3d7e
for n in keys(thicknesses_mm)
	all_data[String(n)].total_thickness = 1e-3 * getfield(thicknesses_mm,n)
end

# ╔═╡ 56a5f3c9-6a41-4326-83ab-2c19d65b3ed0
if is_selected
		for  ki in keys(selected_variables_multi)
			name = string(ki)
			d = all_data[name]
			d.selected_names = getfield(selected_variables_multi , ki)
			fill_data_connector_data!(d)
		end
		
end;

# ╔═╡ b3c96eae-64a5-4245-85b6-b7b994e03ff7
if is_selected
	selected_variables_multi
	@bind locations_data  confirm(multi_values(all_data))
end

# ╔═╡ 735d3901-2dee-40d9-9e74-bb0d71ddfda4
begin # filling thickness data 
	for (_,d) in all_data
		n = length(d.selected_names)
		d.locations = fill(0.0 , n)
	end
	for k in keys(locations_data)
		(projn , sn) =Tuple(split(String(k) , ":"))
		val = getfield(locations_data,k)
		d_i = all_data[projn]
		fl = d_i.selected_names .== sn
		_v = @view d_i.locations[fl] 
		fill!(_v , 1e-3 * val)
	end
end

# ╔═╡ b1f00c55-5ce9-4a0f-a548-a8a9041d02fc
begin 
	is_data_ready = true
	locations_data
	thicknesses_mm
	selected_variables_multi
	cur_proj
	refit
	
	for (_ , d) in all_data
		
		for fi in fieldnames(DataConnector)
			global is_data_ready = is_data_ready && !isempty(getfield(d , fi)) 
			if !is_data_ready 
				display(String(fi) * " is not ready" )
				break
			end
		end	
		global is_data_ready = is_data_ready && all(Base.Fix2(>=,0.0) , d.locations)
		is_data_ready = is_data_ready && d.total_thickness > 0.0
		is_data_ready || break
	end
	is_data_ready

end

# ╔═╡ fc79d398-03d0-4f75-a777-41e2e27fee5e
if is_data_ready 
	(lmin , lmax) = (1e6 , 0.0)
	for (_,d) in all_data
		global lmin , lmax
		t_data= @view d.data_cutted[:,2:end]
		(lmin_cur , lmax_cur) = extrema(t_data)
		(lmin_cur < lmin) && (lmin = lmin_cur) 
		(lmax_cur > lmax) && (lmax = lmax_cur) 
	end
	lam_T_range = (lmin , lmax)
end

# ╔═╡ 16337bb4-438e-4b86-bdc5-b88ef210a960
begin 
	C_poly = IHT.ScaledPolynomial(IHT.BernsteinSymPoly(ntuple(_->density*1000.0 , basis_degree)), xmin = lam_T_range[1], xmax =lam_T_range[2])
	c_x = Cp_Dict[material_name].x 
	c_y = density.*Cp_Dict[material_name].y
	flc = @. (c_x <= lam_T_range[2]) & (c_x >= lam_T_range[1])
	c_x = c_x[flc]
	c_y = c_y[flc]
	if is_optimize_c
		C = IHT.OptimizableVariable(C_poly, flag = is_optimize_c , lb = density*lower_cp_limit , ub = density * upper_cp_limit )
	else

		#C = PW.polyfit!(C_poly, c_x,c_y)#linear_interpolation(c_x ,c_y , extrapolation_bc=Line())
		C =linear_interpolation(c_x ,c_y , extrapolation_bc=Line())
	end
end;

# ╔═╡ c8861fa2-9cb3-4349-9ece-eba285c24eab
begin 
	λ_poly = IHT.ScaledPolynomial(IHT.BernsteinSymPoly(ntuple(_->1.0 , basis_degree)), xmin = lam_T_range[1], xmax =lam_T_range[2])
	
	if is_optimize_lambda
		
		λ = IHT.OptimizableVariable(λ_poly, flag = is_optimize_lambda , lb = lower_lam_limit , ub = upper_lam_limit )
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

# ╔═╡ e06ac758-c5ae-4cf6-84b1-11d51a56f200
initial_distribution = lam_T_range[1]

# ╔═╡ d71fd6d1-167d-40fb-a253-b0982e19c0d0
@bind Ttest Slider(range(lam_T_range[1],lam_T_range[2],100), show_value = true)

# ╔═╡ 1318f293-f606-4a9f-8896-853b87c0665d
begin
	xpoints_number = dp_xpoints
	tpoints_number = dp_tpoints
	#thickness = 1e-3 * thickness_mm
	if is_data_ready
		probls = []
		for (_,d) in all_data
			therm_locations = d.locations
			thickness = d.total_thickness
			inds = sortperm(therm_locations)
			cov = covariance_type(cargs...)
			t = d.data_cutted[: , 1]
			T = d.data_cutted[: , 2:end]
			inv_probl = IHT.SingleInverseProblem(t, T[:,inds], initial_distribution, therm_locations[inds], C,λ, dλdT, thickness, xpoints_number, tpoints_number , covariance=cov , regularization = reg_type())
			push!(probls , inv_probl)
		end
		parallel_probls =IHT.ParallelInverseProblems(Tuple(probls)...)
	end
end;

# ╔═╡ fe228354-5f3f-433d-9f8b-7d888e9c5ecb
begin 
	IHT.set_regularization_multiplier!(parallel_probls , reg_multiplier)
end

# ╔═╡ d33fc044-c487-449f-8f7c-0bcfee603358
parallel_probls

# ╔═╡ 6cd83678-860c-41c8-b3cd-007725f9e01d
begin 
	
	bc_plot = Plots.plot(;plot_common_args...)

	for inv_probl in parallel_probls.problems
		t_dir = OHT.trange(inv_probl.direct_problem)
		Plots.plot!(bc_plot, t_dir ,inv_probl.direct_problem.bc_dwn.(t_dir), label = "lower bc"; plot_common_args...)
		
		Plots.plot!(bc_plot, t_dir ,inv_probl.direct_problem.bc_up.(t_dir), label = "upper bc"; plot_common_args...)
		for (i,c) in enumerate(eachcol(inv_probl.Tdata_measured))
			Plots.plot!(bc_plot, t_dir, c, label="T$(i)" )
		end
		xlabel!(bc_plot,"Time,s")
		ylabel!(bc_plot,"Temperature,s")
	end
end

# ╔═╡ 544201c6-90b3-48dd-9651-69ca3c5c4979
bc_plot

# ╔═╡ dd8fb4fa-fa67-4e26-988d-3ede79bc9540
if use_l_curve
	#res_l = is_after_fit ? res2 : res 
	
	alpha_range = collect(range(extrema(l_curve_alpha)... , length=l_curve_points_number))
	
	alphas = is_l_curve_logscale ? 10.0 .^alpha_range : alpha_range

	#prob_fun = Main.regularization_scan(alphas , Optimization.remake)
	lcurve_probs = IHT.ParallelInverseProblems(Tuple(probls)...)
	
	(_start, _lb, _ub)  =  IHT.fill_starting_vectors(lcurve_probs)
	#@. _start = res_l.u
	
	
	sols = Dict{Int , Any}()
	Threads.@sync for i in 1:length(alphas)
		 Threads.@spawn begin
			p_i = deepcopy(lcurve_probs)
			IHT.set_regularization_multiplier!(p_i , alphas[i])
			optp_i = OptimizationProblem(IHT.discrepancy!, 
										 _start, 
										 p_i,
										 lb = _lb, ub = _ub)
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
	Plots.scatter!(l_curve_plot , l_curve_cov_loss , l_curve_reg_loss ; scale_args..., plot_common_args...)
	xlabel!(l_curve_plot , "Model prediction loss")
	ylabel!(l_curve_plot , "Regularization loss"  )
	
end;

# ╔═╡ e5569bf4-b14f-4806-bbce-38096e2f837b
use_l_curve && reg_plot

# ╔═╡ 985eb56b-c964-4be8-973a-cec2c99df494
use_l_curve && l_curve_plot

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

# ╔═╡ e9e9e16d-0b2c-45ce-aa5b-cdcda6b143f1
begin 
	reg_multiplier
	is_data_ready
	
	if is_fit_on 
	(start, lb, ub)  =  IHT.fill_starting_vectors(parallel_probls)
	optp = OptimizationProblem(IHT.discrepancy!, start, parallel_probls, lb = lb, ub = ub , maxiters=pso_iters )
	res = solve(optp, optimizer())
		IHT.discrepancy!(res.u, parallel_probls)
		if is_after_fit
			new_start = res.u
			disc_before = IHT.discrepancy!(new_start, parallel_probls)
			optp2= OptimizationProblem(IHT.discrepancy!, new_start, parallel_probls , lb = lb, ub = ub , maxiters=pso_iters)
			res2 = solve(optp2, NLopt.LN_COBYLA())
			disc_after = IHT.discrepancy!(res2.u, parallel_probls)
			disc_before - disc_after
		end
	end
	refresh_graph = false
end

# ╔═╡ 4c9a5a5f-5839-4211-b561-7539dfa74a7a
# plotting covariance 
begin 
	refresh_graph
	
	IHT.loss_distribution(parallel_probls)
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
DataFrames.DataFrame(loss_table)

# ╔═╡ f2943cf8-ffb9-4065-a272-1f344488dd0f
begin 
	refresh_graph

	# plotting fitted values 
	(Tmin,Tmax) = (lam_T_range[1],lam_T_range[2])
	Tplot = range(Tmin, Tmax, 100)
	p_fit_lam = Plots.plot(Tplot, λ.(Tplot), label = "inverse", linewidth = 2,linecolor = :green)
	if haskey(L_Dict, material_name) && show_passport
		d_cur = L_Dict[material_name]
		f_range =  (d_cur.x .>= Tmin) .& (d_cur.x .<=Tmax)
		Plots.plot!(p_fit_lam, d_cur.x[f_range],d_cur.y[f_range], label = "passport",marker = :circle;  plot_common_args...)
	end
	xlabel!(p_fit_lam, "Temperature, ᵒC")
	ylabel!(p_fit_lam, "Thermal conductivity, W/(m*K)")
	title!(p_fit_lam, "Fitted VS measured")

	# plotting temeprature distribution 
	p_residual = Plots.plot(;plot_common_args...)
	p_distr = Plots.plot(;plot_common_args...)
	for inv_probl in parallel_probls.problems
		discr = IHT.evaluate_loss(inv_probl)
		Plots.plot!(p_distr , inv_probl.Tdata_evaluated, label = "fitted")
		Plots.plot!(p_distr, inv_probl.Tdata_measured, linestyle = :dash, label = "measured"; plot_common_args...)
		title!(p_distr,"discrepancy = $(discr)")
		xlabel!(p_distr, "Timestep")
		ylabel!(p_distr, "Temperature, ᵒC")
		
		
		# ploting residuals
	
		Plots.plot!(p_residual , inv_probl.residual, label = nothing ; plot_common_args...)
		title!(p_residual, "Residuals")
		xlabel!(p_residual, "Timestep")
		ylabel!(p_residual, "ΔT, ᵒC")
	end
end ;

# ╔═╡ a660bf26-5b50-4910-b0ab-8e453623dc1a
p_residual

# ╔═╡ 538487b4-b2fc-42c2-ba69-663b2ca5b768
begin 
	c_plot = Plots.plot(Tplot,C.(Tplot) ./density, label = "inverse", xlabel = "Temperature",linecolor = :green , ylabel = "Heat capacity, J/K"; plot_common_args...)
	Plots.plot!(c_plot ,c_x ,c_y ./density , label = "passport",marker = :circle , markersize  = 8; plot_common_args... )
end

# ╔═╡ 67763d8f-e1ac-4b1f-978f-4734af1e03ba
ylims!(p_fit_lam, lam_y_scale_region)

# ╔═╡ 042ea29b-63dd-43d0-a20f-68807c5f7cd4
p_distr

# ╔═╡ 27031733-7ade-4bda-b464-5a9ade0b2950
begin 
	name_poly =  material_name
	x_fit = L_Dict[name_poly].x
	www =   (d_cur.x .>= Tmin) .& (d_cur.x .<=Tmax)
	x_fit = x_fit[www]
	y_fit = L_Dict[name_poly].y[www]
	(k,_,rrr) =  PW.scale_x_to_ξ(x_fit)
	pp = PW.polyfit(PW.BernsteinSymPoly{basis_degree, Float64} , k , y_fit)
	ppp = Plots.plot(x_fit,pp.(k), label = nothing)
	Plots.plot!(ppp,x_fit,y_fit, marker=:circle, label = nothing; plot_common_args...)
	Plots.scatter!([Ttest] ,[lam_test] )
	title!(ppp,"""Polynomial fitting of passport,
		   r²=  $(norm(y_fit .- pp.(k))/norm(y_fit))
		   """)
end

# ╔═╡ 4cf77d64-a720-41da-a9c2-5d875ea00135
begin 
	is_write_file
	#save_folder = select_folder_gtk()
	Tsave = range(Tmin, Tmax, 100)
	mat = hcat(Tsave, λ.(Tsave))
	fullfile_name = joinpath(save_folder,default_save_file)
	CSV.write(fullfile_name,Tables.table(mat, header = ["T", "lambda"] ), delim = " ")
end

# ╔═╡ a2a84394-1f18-48e3-a4dc-c03f64a3a1c0
begin
	refresh_graph
	p_hist = Plots.histogram()
	#Plots.histogram(inv_probl.residual[:,1] )
	for inv_probl in parallel_probls.problems
	for (i,c) in enumerate(eachcol(inv_probl.residual))
		bb = range(minimum(c),maximum(c),20)
		Plots.histogram!(p_hist,c,bins = bb, alpha=0.5, label="T$(i)"; plot_common_args...)
	end
	end
	p_hist
end

# ╔═╡ a1c795e4-1be0-46e0-b6a2-4eda990fd65e
function time_range_selector(all_data::AbstractDict)
	N = length(all_data)
	PlutoUI.combine() do Child
		@htl("""
		<h6>Time range, tmin - tmax , s</h6>
		<ul>
		$([
			@htl("<li>$(k): $(Child( k, RangeSlider( minimum(d.data[:,1]) : maximum(d.data[:,1]) )))</li>")
			
			for (k , d) in all_data
		])
		</ul>
		""")

	end
end

# ╔═╡ df99b7b7-5a0d-4b71-bf5c-415b5cfa3b2d
if is_selected
	selected_variables_multi
	@bind time_region  confirm(time_range_selector(all_data))
end

# ╔═╡ 5be15388-cea2-4884-b8de-bff5be64e506
if is_selected 
	raw_data_plot = Plots.plot(;plot_common_args...)
		
	for (k,d) in all_data
		rng = getfield(time_region , Symbol(k))
		(tmin , tmax) = extrema(rng)
		cut_data!(d , tmin , tmax)
		!isempty(d.data_cutted) || continue
		leg_str = first(eachsplit(k , "_"))
		for (i , c) in enumerate(eachcol(d.data_cutted)[2:end])
			Plots.plot!(raw_data_plot , d.data_cutted[:,1] , c , label = leg_str * ":" * d.selected_names[i]  , legend_position = :best   , linestyle = :auto;plot_common_args...)
		end
	end
	raw_data_plot
end

# ╔═╡ dbd1d788-120b-4678-bf3b-7541b9ea7341
begin
	time_region
	let d = all_data[show_data_table]
		tmin_c = t_cutted_min(d)
		tmax_c = t_cutted_max(d)
		tbl = WP.to_table(d.project ; names = d.selected_names , tmin = tmin_c , tmax = tmax_c)
	end
end

# ╔═╡ 87ae11cd-8c4f-4835-9a9f-852447a45c97
#=write(joinpath(raw"D:\JULIA\JULIA_DEPOT\dev\InverseHeatTransfer.jl\test\test_data\binary_files","dasfsdf.json"), JSON2.write(Dict(:T1=>2.0, :T2=>3.0))) =#

# ╔═╡ 58eb27e1-38c5-4a90-9e80-61f6213aa721
#= dd = JSON2.read(read(joinpath(raw"D:\JULIA\JULIA_DEPOT\dev\InverseHeatTransfer.jl\test\test_data\binary_files","dasfsdf.json"), String)) =#

# ╔═╡ 304091a8-4206-49c2-9c98-cffb18a0e906
# ╠═╡ disabled = true
#=╠═╡
begin
	all_data =Dict{String , DataConnector}()
	for n in  cur_proj
		selected_proj = projects[n]
		all_data[selected_proj.name] = DataConnector(collect(keys(selected_proj.data)) , String[] , Matrix{Float64}(undef,0,0) , Matrix{Float64}(undef,0,0) , Float64[] , 0.0, selected_proj)
	end
end;
  ╠═╡ =#

# ╔═╡ 4c41953d-1625-4383-9e57-545ab7f4c0e5
all_data = DC.DataSelectorsGroup(projects)

# ╔═╡ Cell order:
# ╠═a17fe1fe-5542-454b-b45e-942ac52b6f1a
# ╠═2bb45200-d7ba-47a3-b82d-9c147b5d7601
# ╠═80479ac3-897a-4aa7-8ce9-977fa4912bc6
# ╠═6d2a48d4-e458-4edd-b22a-7b7dae9f492c
# ╠═a6d591eb-20f2-476f-af19-9b0084ddf929
# ╠═35958e8a-eb7a-4eff-89a0-f9c04aff2a37
# ╠═438a3909-9367-4660-a3ca-bd1786ab6016
# ╟─db671921-13dc-497b-81e5-dcb4da0695f9
# ╠═450fb200-eec6-4e96-9ebd-81453c015830
# ╠═a407a99b-b40c-436c-a2a0-af2e19b347b8
# ╟─6e062bd9-d20c-4e1d-b772-328bec8859ea
# ╠═81c2f93b-05b1-4eb0-9919-4ef76ecad233
# ╟─17fbc55f-12a8-431e-ac00-b28304f2eb6c
# ╠═4c41953d-1625-4383-9e57-545ab7f4c0e5
# ╠═0a4ce046-5b43-425d-a3f2-da8d9867b181
# ╠═dcf0034b-e405-4f27-b853-acb7a58b9cc6
# ╟─4d27d75e-5e95-43c9-9e77-b9ae3d6c4bb9
# ╟─28b91e49-996c-4abc-9db4-30b515444aab
# ╠═054d932d-12be-4538-ab34-b5d3f465bf0f
# ╠═304091a8-4206-49c2-9c98-cffb18a0e906
# ╠═2e676f8c-909a-4fb7-b35e-57d08f802df7
# ╠═0177413c-0f89-4935-9963-5aeebd333b9a
# ╠═c30f95e5-93eb-4ab7-bc84-4bfe7092cf47
# ╟─56a5f3c9-6a41-4326-83ab-2c19d65b3ed0
# ╟─bc6ecff4-2ade-4436-9630-be573eb1ea04
# ╟─df99b7b7-5a0d-4b71-bf5c-415b5cfa3b2d
# ╟─5be15388-cea2-4884-b8de-bff5be64e506
# ╟─b3c96eae-64a5-4245-85b6-b7b994e03ff7
# ╟─6c4cb363-3bda-4dfa-8499-8201a013895d
# ╠═dbd1d788-120b-4678-bf3b-7541b9ea7341
# ╟─735d3901-2dee-40d9-9e74-bb0d71ddfda4
# ╟─ca69ae4b-8428-4643-93a4-c4c1687b3d7e
# ╠═fd47bc58-8a0c-4dd7-875c-bcc80a21e64e
# ╠═b1f00c55-5ce9-4a0f-a548-a8a9041d02fc
# ╟─62069546-44fb-4a77-986d-f03624719e29
# ╟─fc5c1209-26ec-41d7-9238-e57d18330de1
# ╟─2b3a1a65-8d3c-424e-a6cf-0e96646795f4
# ╟─8f65e989-1abe-4689-9c3f-fdb6cabd7eae
# ╟─a46970c8-7c91-4795-83a9-41c56a7ca399
# ╟─5843fcfb-4e0e-480d-b263-e3f3ad6a7ac3
# ╟─7d37019a-34bf-43ca-aa81-8b6c29191473
# ╟─4f3e1899-541d-4809-810c-f2e6b7ca1ad1
# ╟─2bb61e43-384c-48ce-8a55-843726bf3f05
# ╟─baffb68a-fb52-4bac-b484-4a215108aaed
# ╟─16337bb4-438e-4b86-bdc5-b88ef210a960
# ╟─4ca95124-8a7a-4e4f-9e65-ff7b2adf35a5
# ╟─c8dc4f9d-a549-4dc2-82bd-38ffe949ea55
# ╟─bfa23359-8bac-4db0-bac1-0885ebe8ec4b
# ╠═a660bf26-5b50-4910-b0ab-8e453623dc1a
# ╠═538487b4-b2fc-42c2-ba69-663b2ca5b768
# ╟─67763d8f-e1ac-4b1f-978f-4734af1e03ba
# ╟─68ed85b2-7308-4ea7-b696-0f1951219592
# ╟─9b35c13b-24ba-4c43-a621-a3f8ba45fe4a
# ╟─33473b7a-e22f-4333-a2f2-374778c0d603
# ╟─672119a2-7a47-4813-a2d3-e0c15ee63491
# ╟─c8cd1797-b631-4fdf-a51f-67e22c86c55a
# ╟─eb70f6d9-301f-4e0f-91a2-0e911b1a6dd9
# ╟─10f7c083-a697-4b7e-87ea-de2791ed1a30
# ╟─3f8c5c6a-f902-4b1e-abc1-7b572c8d2505
# ╟─e5569bf4-b14f-4806-bbce-38096e2f837b
# ╟─985eb56b-c964-4be8-973a-cec2c99df494
# ╟─0a0324df-430f-4ac0-857b-4da6a7dca138
# ╟─fe228354-5f3f-433d-9f8b-7d888e9c5ecb
# ╟─f3e24d8a-e35d-41de-92e6-0df290d3503c
# ╟─4c9a5a5f-5839-4211-b561-7539dfa74a7a
# ╟─95dbc55f-ab5a-4828-a1e2-9a0c9a9ec19b
# ╟─fa72774e-040a-4bc3-a759-5eb68c243fb4
# ╟─3281dd3c-0453-4098-a6ed-4d771717033b
# ╟─bbc2d0df-28bf-4754-bbdf-bece0bbfe76b
# ╟─2a1ca349-3732-443e-bc48-5be611a5d91f
# ╟─d02ec3fd-1b0d-4bc3-92e9-adedc8bf2c8c
# ╟─fc79d398-03d0-4f75-a777-41e2e27fee5e
# ╟─c8861fa2-9cb3-4349-9ece-eba285c24eab
# ╟─e06ac758-c5ae-4cf6-84b1-11d51a56f200
# ╟─1318f293-f606-4a9f-8896-853b87c0665d
# ╟─d33fc044-c487-449f-8f7c-0bcfee603358
# ╟─6cd83678-860c-41c8-b3cd-007725f9e01d
# ╟─544201c6-90b3-48dd-9651-69ca3c5c4979
# ╟─f2943cf8-ffb9-4065-a272-1f344488dd0f
# ╟─042ea29b-63dd-43d0-a20f-68807c5f7cd4
# ╟─a2a84394-1f18-48e3-a4dc-c03f64a3a1c0
# ╟─9cd8c4c8-7d11-4fb0-92d5-439702aa9496
# ╠═e9e9e16d-0b2c-45ce-aa5b-cdcda6b143f1
# ╟─aea6fc8f-5ee2-40a3-aebd-942c45eec0d6
# ╠═a8f5cc07-b2bd-445c-8797-ffe78a25641b
# ╠═bce5f7ae-49c5-4520-a0f6-6215e5078674
# ╠═dd8fb4fa-fa67-4e26-988d-3ede79bc9540
# ╠═cf8665e2-07d9-4464-b833-0531205408e9
# ╠═2614e5ad-0d45-4046-abfd-053dd872c17c
# ╠═1c32440c-32af-414c-9a7a-d6cc83685f7c
# ╠═84e23b57-e5a4-4bf0-97ef-6b41b502b4e6
# ╠═c2564f33-dcf5-45fe-a462-70fe105bcf0a
# ╟─27031733-7ade-4bda-b464-5a9ade0b2950
# ╟─f5b13ec5-98f9-48bb-ad95-7dbd87e11f7b
# ╟─d71fd6d1-167d-40fb-a253-b0982e19c0d0
# ╟─fb2937d0-738b-4329-aef5-f3af1d17497a
# ╟─efa9120f-5a45-41a0-9132-dc26f967fec3
# ╟─023dc25f-6cf8-4802-83b4-77d1827cd2a7
# ╟─23902c7d-5973-4868-8488-e0c7634573c4
# ╟─6800ae42-c5b1-4d1e-84fb-a03015bf138f
# ╟─ea412a80-0bf7-4ac6-9b90-8530a4c26008
# ╟─4cf77d64-a720-41da-a9c2-5d875ea00135
# ╠═2bf91a7e-0da8-48e8-a779-962be2e7c03b
# ╟─fba5bc8b-25cb-406b-aa13-f06c591e08c9
# ╟─5135015e-0fcf-484a-86ad-3cb7e40849bd
# ╟─fd51a6c8-6569-4bfd-86ac-883c648fe6d9
# ╠═6cd4f554-7242-4e97-b4cb-8549e3b70139
# ╟─6d0d7cb4-45fc-405f-8a5d-cf10ad5e380b
# ╟─646aebc7-b3db-443f-888f-800d444b4fa3
# ╟─60458134-1de0-475e-ba64-24b6f33a2980
# ╟─96b4c7c4-30a8-407c-9bd4-245f0ee0d9b2
# ╟─a1c795e4-1be0-46e0-b6a2-4eda990fd65e
# ╠═87ae11cd-8c4f-4835-9a9f-852447a45c97
# ╠═58eb27e1-38c5-4a90-9e80-61f6213aa721
