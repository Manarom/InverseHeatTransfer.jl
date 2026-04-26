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
	using StaticArrays , Plots  , Revise , PrettyTables
	using PlutoUI ,  OrderedCollections
	using HypertextLiteral
	import InverseHeatTransfer
	using Dates
end

# ╔═╡ 5807712b-5d26-49c8-ab65-dac167ebad7b
begin 
	projects_save_path_ref =  Ref("") # projects saving path
	projects_save_name_ref =  Ref("") # projects saving name

	data_selection_save_path_ref = Ref("") # data selection saving path
	data_selection_save_name_ref = Ref("") # data selection saving name
	
	default_data_fodler = joinpath(@__DIR__, "..","src","data_utils","binary_files")
	source_path = joinpath(@__DIR__,"..","src")

	includet(joinpath(@__DIR__ , "CustomPlutoFunctions.jl"))
end;

# ╔═╡ 35958e8a-eb7a-4eff-89a0-f9c04aff2a37
begin 
	IHT = InverseHeatTransfer
	PW = IHT.ScaledPolynomials
	OHT = IHT.OneDHeatTransfer
	WP = IHT.WinPos
	DC = IHT.DataConnector
	PF = Main.CustomPlutoFunctions
	plot_common_args = (grid = true, gridlinewidth=3, gridstyle = :dot,minorgrid=true, box = :on, linewidth = 3 , xlabel = "time , s" , ylabel = "temperature , ᵒC");
end;

# ╔═╡ db671921-13dc-497b-81e5-dcb4da0695f9
md""" ## Loading experimental data from winpos project"""

# ╔═╡ 450fb200-eec6-4e96-9ebd-81453c015830
md" Load data from : $(@bind input_data_type Select([:winpos, :ascii , :hdf5_winpos , :hdf5_data_selector] , default = :winpos))"

# ╔═╡ 6e062bd9-d20c-4e1d-b772-328bec8859ea
begin 
	md""" #### working folder $(@bind working_folder TextField(90, default = realpath(default_data_fodler))) """
end

# ╔═╡ 7370790d-cc53-4159-84b1-20e8839a8fc8
if @isdefined(working_folder) && isdir(working_folder) 
	all_hdf5_files = [d for d in readdir(working_folder) if contains(d , ".hdf5") ]
	if isempty(all_hdf5_files)
		md" **There is no hdf5 files in the folder**"
	else
		md" **Select data selection project** $(@bind data_selection_name Select(all_hdf5_files))"

	end
else
	@isdefined(working_folder) ? md"**Not a folder:**" : md"**Not a folder: $(working_dir)**" 
end

# ╔═╡ 81c2f93b-05b1-4eb0-9919-4ef76ecad233
begin 
	projects = nothing 
	try 
		global projects =
			if input_data_type == :winpos
				WP.load_from_winpos_folder(working_folder)
			elseif input_data_type == :ascii
				WP.load_from_ascii_folder(working_folder , name_matcher = "T")
			elseif input_data_type == :hdf5_winpos
				_full_name = joinpath(working_folder , data_selection_name)
				WP.load_from_hdf5(_full_name , WP.WinPosProjectsGroup)
			elseif input_data_type == :hdf5_data_selector
				_full_name = joinpath(working_folder , data_selection_name)
				_ds = WP.load_from_hdf5(_full_name , DC.DataSelectorsGroup)
				convert(WP.WinPosProjectsGroup, _ds)
			else
				error("Unknown data format")
			end
	catch er 
		md" Unable to load file due to $(er)"
	end
	is_any_projects = !isnothing(projects)	
end;

# ╔═╡ 17fbc55f-12a8-431e-ac00-b28304f2eb6c
is_any_projects && md""" #### select projects: 
$(@bind selected_projects confirm(MultiSelect(collect(keys(projects)) , default = ["10ks" , "2ks"]))) 
"""

# ╔═╡ d9bac527-2858-4536-b35a-ecd03fb11ec8
 md"""
	**$((is_projects_selected = (is_any_projects && !isempty(selected_projects))) ? "Several projects selected " : "Select projects to proceed" )**
	"""

# ╔═╡ 46dc9aed-cbbe-4b8c-a6e3-9a5207bee10c
md"""
	#### **Save selected projects group to file ? :** $(@bind save_selected_projects_trigger PlutoUI.CheckBox(false))

	"""

# ╔═╡ ad3152ec-d481-4b65-856b-f7c1987d379e
@bind project_saving_type Select(["hdf5" , "winpos folder"] , default = "hdf5")

# ╔═╡ f9bf69a8-9854-4c46-8c1c-e4af8ed176d8
is_any_projects && @bind projects_save_name_path PF.multi_values(PlutoUI.TextField , ("path" , "name") , title = "Path/Name" , defaults=(working_folder , "projects"))

# ╔═╡ c9fa68be-e2c4-4c4c-a6cc-f9f064a0a4c8
if is_any_projects 
	projects_save_path_ref[] = projects_save_name_path.path 
	projects_save_name_ref[] = projects_save_name_path.name
end;

# ╔═╡ 6e7839ca-da24-4a3e-927f-b035729a4cb7
begin  
	
	if save_selected_projects_trigger && is_projects_selected
		
		isdir(projects_save_path_ref[]) || mkdir(projects_save_path_ref[])
		is_hdf5 = project_saving_type == "hdf5"
	
		selected_projects_group = WP.WinPosProjectsGroup(projects.name, OrderedDict([p for p in projects if first(p) ∈  selected_projects]...) , working_folder)
		_ffile	= joinpath(projects_save_path_ref[], projects_save_name_ref[])
		 is_hdf5 &&  (_ffile *= ".hdf5") 
		if is_hdf5 
			WP.export_to_hdf5(selected_projects_group , _ffile)
		else
			isdir(_ffile) || mkdir(_ffile)
			WP.write_to_winpos_folder(selected_projects_group , _ffile)
		end
		"✅ Saved to $(project_saving_type)-file $(_ffile) at $(Dates.format(now(), "HH:MM:SS"))"
	else
		"😖 Not saved at $(Dates.format(now(), "HH:MM:SS"))"
	end
end

# ╔═╡ 3aa5a908-c4eb-4f6a-899e-c7ba2d29fa01
if is_any_projects && is_projects_selected 
	all_data = DC.DataSelectorsGroup([ p for p in projects if first(p) ∈  selected_projects]...)
end;

# ╔═╡ c30f95e5-93eb-4ab7-bc84-4bfe7092cf47
is_projects_selected && @bind selected_variables_multi  confirm(PF.multi_values_table(PlutoUI.MultiSelect , all_data.d , WP.all_names))

# ╔═╡ 3126537e-cd1d-476d-9f0d-84b37ac15678
is_sensors_selected = @isdefined(selected_variables_multi) && !any(isempty , selected_variables_multi);

# ╔═╡ 56a5f3c9-6a41-4326-83ab-2c19d65b3ed0
if is_sensors_selected
	DC.unselect!(all_data)
	foreach(pairs(selected_variables_multi)) do (ki , ni)
		DC.select!(all_data , String(ki) , ni)
	end
	DC.fill_data_for_selected!(all_data)
end;

# ╔═╡ bc6ecff4-2ade-4436-9630-be573eb1ea04
if is_sensors_selected
	@bind thicknesses_mm confirm(PF.multi_values_table(PlutoUI.NumberField , collect(keys(all_data)) , default = 0.0 , fontsize = 20 ,  column_names = ("Project" , " Thickness ") ,  title = "Thickness , mm"))
end

# ╔═╡ 999cefa4-3513-4c4f-86f1-49d4d55c66f8
if is_sensors_selected
	thicknesses_mm
	is_thickness_set = true
	try 
		foreach(pairs(thicknesses_mm)) do (k,v)
			d_i = all_data[String(k)]
			DC.set_thickness!(d_i , 1e-3*v)
		end
		global is_thickness_set = true
	catch er 
		md" **$(er)**"
		global is_thickness_set = false
	end
end;

# ╔═╡ c6764c6c-3504-42f8-9140-d40e7d050000
is_projects_selected && @bind show_raw_data_table Select(collect(keys(all_data)))

# ╔═╡ 1730c035-2c46-4267-86ee-0cfc60073f96
md" Show table $(@bind is_show_selected_data_table CheckBox(false))"

# ╔═╡ 3053a02a-15a2-444f-b555-a44876bdf0a6
if is_projects_selected && is_sensors_selected
	raw_raw_data_plot = Plots.plot(;plot_common_args...)
	try
	_data_combined_raw = DC.selected_data(all_data[show_raw_data_table])
	raw_raw_headers = ["t" , DC.selected_names(all_data[show_raw_data_table])...] 

	
	raw_raw_t = @view _data_combined_raw[:,1] 
	raw_raw_labels = @view raw_raw_headers[2:end]
	for (i,c) in enumerate(eachcol(@view _data_combined_raw[:,2:end]))

		Plots.plot!(raw_raw_data_plot , raw_raw_t ,  c; label = raw_raw_labels[i], plot_common_args... )
	end
			catch er 
		md" $(er)"
	end
end

# ╔═╡ c1e54cfc-da70-4e29-8da7-19227cd56e6d
if is_show_selected_data_table && is_projects_selected && is_sensors_selected

	
	@htl("""
	<div style="max-height: 200px; overflow-y: auto; border: 1px solid #ccc;">
	    $(pretty_table(HTML , _data_combined_raw , column_labels = raw_raw_headers, top_left_string ="Temeratures for $(show_raw_data_table)"))
	</div>
	""")

end

# ╔═╡ b8d5a2ee-5a5e-462f-9588-c7d910f446e8
is_projects_selected && raw_raw_data_plot

# ╔═╡ 6c4cb363-3bda-4dfa-8499-8201a013895d
is_projects_selected && @bind show_data_table Select(collect(keys(all_data)))

# ╔═╡ 29a0e7ad-c939-4dbc-a724-6afbad32ff5f
md" Show table $(@bind is_show_all_data_table CheckBox(false))"

# ╔═╡ c3f0a382-f886-41ae-86ae-e5dbebc16214


# ╔═╡ d0118c80-2ece-4185-bd8c-76983b9954f2
"""$(Dates.format(now(), "HH:MM:SS dd:mm:yyyy"))"""

# ╔═╡ 6da6a4dd-2bf4-463c-91a0-5e73b3e3a1ef
@bind comment_to_sample html"""
<textarea 
	rows='10' 
	cols='50' 
	placeholder=
	'Enter sample description here... 
					\n $(Dates.format(now(), "HH:MM:SS"))'
	style='font-size: 20px; font-family: "Courier New", monospace; line-height: 1.5; padding: 10px;'
></textarea>
"""

# ╔═╡ 77e39e2f-4a92-4add-944e-b8f7d4d139df
comment_to_sample

# ╔═╡ c9717703-5218-451d-9a80-a4ddb79b929d

@bind data_selection_save_pathname PF.multi_values(PlutoUI.TextField, ("path" , "name") , title = "Data selection path/name" , defaults=(working_folder , "data_selection"))

# ╔═╡ 8e0fff6a-4a02-4e6b-a017-477a46269918
md""" Resave data selection $(@bind resave_selection_trigger Button("resave"))"""

# ╔═╡ 85c827cd-c4c8-4c0a-b790-86b2ea070e1d
begin 
	data_selection_save_path_ref[] = data_selection_save_pathname.path
	data_selection_save_name_ref[] = data_selection_save_pathname.name
end;

# ╔═╡ 12699165-eb53-4bca-b177-8326a0a72aed
md"""

**Save data selection ?** $(@bind data_selection_save_trigger PlutoUI.CheckBox(false))
"""

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

# ╔═╡ c40ec284-b81a-4039-bb33-238de0ca09e4
function thermocouples_locations(d::DC.DataSelectorsGroup)
		PF.multi_values_table(NumberField , DC.joined_selected_names(d) , title = "Thermocouple locations, mm" , default= 0.0)
end

# ╔═╡ b3c96eae-64a5-4245-85b6-b7b994e03ff7
if is_sensors_selected && is_thickness_set
	selected_variables_multi 
	
	@bind locations_data  confirm(thermocouples_locations(all_data))
end

# ╔═╡ ea647fc8-37c7-4726-980d-42f97ffc02e3
if is_sensors_selected && is_thickness_set
	
	locations_data

	try 
		foreach(pairs(locations_data)) do (k,v)
			(name , couple ) = String.(split(String(k) , ":"))
			d_i = all_data[name]
			DC.set_location!(d_i , couple , 1e-3*v)
		end
	catch ex
		md" $(ex.msg)"
	end
end

# ╔═╡ dd0ecd31-3cf5-4d4d-b3f8-125b5f899c29
function range_constructor(d::DC.DataSelector; npoints::Int = 3000) 
	(tmin, tmax) = DC.default_tmin_tmax(d)
	return range(round(tmin), round(tmax) , npoints)
end

# ╔═╡ 6b4d07af-a40b-4ed1-a610-2a433311c94e
function time_range_selector(all_data::DC.DataSelectorsGroup)

	PF.multi_values_table(PF.PlutoUI.RangeSlider , all_data.d , range_constructor ; title = "Select time range" , show_value = true)
end

# ╔═╡ df99b7b7-5a0d-4b71-bf5c-415b5cfa3b2d
if is_sensors_selected && is_thickness_set
	selected_variables_multi
	@bind time_region  time_range_selector(all_data)
end

# ╔═╡ b31a7c52-7c4b-480f-84d8-fcb1c71dca1b
if is_sensors_selected && is_thickness_set 
	locations_data
	time_region
	thicknesses_mm
	
	selected_plot = Plots.plot(;plot_common_args...)
	_data_combined = DC.combine_selected_data(all_data[show_data_table])
	_t = _data_combined.time_data
	for (i,c) in enumerate(eachcol(_data_combined.temperatures))
		Plots.plot!(selected_plot , _t , c ; label=_data_combined.selected_names[i]  , plot_common_args...)
	end
	title!(selected_plot , show_data_table)
	selected_plot
end

# ╔═╡ c69d07a3-ba0a-48a2-b296-1944b8cd322c
if is_show_all_data_table && is_sensors_selected && is_thickness_set 
	@htl("""
<div style="max-height: 300px; overflow-y: auto; border: 1px solid #ccc;">
    $(pretty_table(HTML , hcat(_data_combined.time_data, _data_combined.temperatures) , column_labels = ["t" , _data_combined.selected_names...] , top_left_string ="Temeratures for $(show_data_table)"))
</div>
""")
end

# ╔═╡ 359994cc-01e8-453d-b3e3-a878ffe466eb
if is_sensors_selected && is_thickness_set
	time_region
	for (k , trange) in pairs(time_region)
		d_i = all_data[String(k)]
		(_tmin , _tmax) = extrema(trange)
		DC.set_time_region!(d_i , tmin = _tmin , tmax = _tmax)
	end
end

# ╔═╡ 35b5c27e-921f-4fe7-90bc-3b327d9150fe
if is_sensors_selected && is_thickness_set 
	
	locations_data
	time_region
	thicknesses_mm
	
	table_data = Matrix{Any}(undef , (length(all_data.d) , 5))
	for (i , (k , d)) in enumerate(all_data.d)
		table_data[i, 1]  = k
		table_data[i, 2] = 1e3*DC.thickness(d)
		table_data[i, 3] = [ Pair(v,k) for (k,v) in zip(1e3*DC.sensors_locations(d) , DC.selected_names(d))]
		table_data[i, 4] = DC.tmin(d)
		table_data[i, 5] = DC.tmax(d)
	end

	pretty_table(HTML , table_data , column_labels = ["name"," h" , "Locs","tmin", "tmax"] , top_left_string ="Sample properties")

end

# ╔═╡ 5be15388-cea2-4884-b8de-bff5be64e506
if is_sensors_selected && is_thickness_set 
	locations_data
	time_region
	thicknesses_mm

	
	raw_data_plot = Plots.plot(;plot_common_args...)
	for (k , d_i)  in DC.selected_data_cutted_with_keys(all_data)
		isempty(d_i.data) && continue
		_t = d_i.data[:,1]
		_names = d_i.names
		CN = size(d_i.data, 2)
		for (i, c) in enumerate(eachcol(d_i.data)[2:CN])
			Plots.plot!(raw_data_plot , _t , c , label ="$(k) : $(_names[i + 1])" ; linestyle = :auto , plot_common_args...)
		end
	end
	raw_data_plot
end

# ╔═╡ b1f00c55-5ce9-4a0f-a548-a8a9041d02fc
begin 
	is_data_ready = true
	try 
		locations_data
		thicknesses_mm
		time_range_selector
		selected_variables_multi
		for (k , d_i) in all_data.d
		 	DC.combine_selected_data(d_i)
		end
		global is_data_ready = true
		md" **Data is ready**"
	catch er
		global is_data_ready = false
		md" Data is not ready due to **$(er)**"
	end
end

# ╔═╡ 7a0abf93-e203-40b2-bcae-3e50c9de5eba
if is_data_ready &&  data_selection_save_trigger
		resave_selection_trigger
		isdir(data_selection_save_path_ref[]) || mkdir(data_selection_save_path_ref[])
	
		_f_f_file	= joinpath(data_selection_save_path_ref[], data_selection_save_name_ref[]*".hdf5")
		WP.export_to_hdf5(all_data , _f_f_file , group_name =data_selection_save_name_ref[] )
	
		"✅ Data selection saved to hdf5-file $(_f_f_file) at $(Dates.format(now(), "HH:MM:SS"))"
end

# ╔═╡ 2390e2fe-1bbe-4a21-a0fb-199cc314a29b
function comment_block(d::DC.DataSelectorsGroup)
	comments_dict["Date"] = """$(Dates.format(now(), "HH:MM:SS dd:mm:yyyy"))"""
	ks = collect(keys(comments_dict))
	def_vals = collect(values(comments_dict))
	PF.multi_values(PlutoUI.TextField , 
				ks , 
				default_values = ntuple(_->50 , length(comments_dict)) , 
				defaults = def_vals)

end

# ╔═╡ b7d63214-f644-42b9-b1c6-60edf7afd903
comment_block(all_data)

# ╔═╡ 5eee07e4-05d6-4125-80ac-fa0777f45dda
#= const comments_dict = OrderedDict(
		"Material"=>DaAn((PUITF,  , "" ) ),
		"Date"=>DaAn( (PUITF , """$(Dates.format(now(), "HH:MM:SS"))""",
		"User"=>"Родин Н.В.",
		"TC type"=>"K",
		"Cement" => "DKS8",
		"Methodics" => "",
		"DataDates" => ""
)=#

# ╔═╡ 77f3b5db-d76e-4518-9354-42d7ebfefa4b
DaAn = NamedTuple{(:type , :default_value , :default)}

# ╔═╡ 0bf33e66-036c-4a6b-aaa7-9b1f60ae5255
PUITF = TextField;

# ╔═╡ Cell order:
# ╟─a17fe1fe-5542-454b-b45e-942ac52b6f1a
# ╟─5807712b-5d26-49c8-ab65-dac167ebad7b
# ╟─35958e8a-eb7a-4eff-89a0-f9c04aff2a37
# ╟─db671921-13dc-497b-81e5-dcb4da0695f9
# ╟─450fb200-eec6-4e96-9ebd-81453c015830
# ╟─6e062bd9-d20c-4e1d-b772-328bec8859ea
# ╟─7370790d-cc53-4159-84b1-20e8839a8fc8
# ╟─81c2f93b-05b1-4eb0-9919-4ef76ecad233
# ╟─17fbc55f-12a8-431e-ac00-b28304f2eb6c
# ╟─d9bac527-2858-4536-b35a-ecd03fb11ec8
# ╟─46dc9aed-cbbe-4b8c-a6e3-9a5207bee10c
# ╟─ad3152ec-d481-4b65-856b-f7c1987d379e
# ╠═f9bf69a8-9854-4c46-8c1c-e4af8ed176d8
# ╟─c9fa68be-e2c4-4c4c-a6cc-f9f064a0a4c8
# ╟─6e7839ca-da24-4a3e-927f-b035729a4cb7
# ╟─3aa5a908-c4eb-4f6a-899e-c7ba2d29fa01
# ╠═c30f95e5-93eb-4ab7-bc84-4bfe7092cf47
# ╟─3126537e-cd1d-476d-9f0d-84b37ac15678
# ╟─56a5f3c9-6a41-4326-83ab-2c19d65b3ed0
# ╟─bc6ecff4-2ade-4436-9630-be573eb1ea04
# ╟─999cefa4-3513-4c4f-86f1-49d4d55c66f8
# ╟─3053a02a-15a2-444f-b555-a44876bdf0a6
# ╟─c1e54cfc-da70-4e29-8da7-19227cd56e6d
# ╟─1730c035-2c46-4267-86ee-0cfc60073f96
# ╟─c6764c6c-3504-42f8-9140-d40e7d050000
# ╟─b8d5a2ee-5a5e-462f-9588-c7d910f446e8
# ╟─b3c96eae-64a5-4245-85b6-b7b994e03ff7
# ╟─ea647fc8-37c7-4726-980d-42f97ffc02e3
# ╟─df99b7b7-5a0d-4b71-bf5c-415b5cfa3b2d
# ╟─6c4cb363-3bda-4dfa-8499-8201a013895d
# ╟─b31a7c52-7c4b-480f-84d8-fcb1c71dca1b
# ╟─29a0e7ad-c939-4dbc-a724-6afbad32ff5f
# ╟─c69d07a3-ba0a-48a2-b296-1944b8cd322c
# ╟─359994cc-01e8-453d-b3e3-a878ffe466eb
# ╟─35b5c27e-921f-4fe7-90bc-3b327d9150fe
# ╟─5be15388-cea2-4884-b8de-bff5be64e506
# ╟─b1f00c55-5ce9-4a0f-a548-a8a9041d02fc
# ╠═c3f0a382-f886-41ae-86ae-e5dbebc16214
# ╠═d0118c80-2ece-4185-bd8c-76983b9954f2
# ╠═6da6a4dd-2bf4-463c-91a0-5e73b3e3a1ef
# ╠═77e39e2f-4a92-4add-944e-b8f7d4d139df
# ╟─c9717703-5218-451d-9a80-a4ddb79b929d
# ╟─8e0fff6a-4a02-4e6b-a017-477a46269918
# ╟─85c827cd-c4c8-4c0a-b790-86b2ea070e1d
# ╟─12699165-eb53-4bca-b177-8326a0a72aed
# ╟─7a0abf93-e203-40b2-bcae-3e50c9de5eba
# ╟─2bf91a7e-0da8-48e8-a779-962be2e7c03b
# ╟─fd51a6c8-6569-4bfd-86ac-883c648fe6d9
# ╟─c40ec284-b81a-4039-bb33-238de0ca09e4
# ╠═6b4d07af-a40b-4ed1-a610-2a433311c94e
# ╠═dd0ecd31-3cf5-4d4d-b3f8-125b5f899c29
# ╠═2390e2fe-1bbe-4a21-a0fb-199cc314a29b
# ╠═b7d63214-f644-42b9-b1c6-60edf7afd903
# ╠═5eee07e4-05d6-4125-80ac-fa0777f45dda
# ╠═77f3b5db-d76e-4518-9354-42d7ebfefa4b
# ╠═0bf33e66-036c-4a6b-aaa7-9b1f60ae5255
