
module DataConnector
	include("WinPos.jl")
	using StaticArrays , OrderedCollections , RecipesBase , Reexport
	using HDF5 , Dates
	@reexport using .WinPos
	export SampleProperties , DataSelector , DataSelectorsGroup ,
	 combine_selected_data 
	const DEFAULT_MEASUREMENTS_SPECIFICATION = OrderedDict(
				"Material"=> "RBSN",
				"Date"=>  """$(Dates.format(now(), "HH:MM:SS dd:mm:yyyy"))""",
				"User"=> "Родин Никита",
				"TCtype"=> "K" ,
				"Cement" =>  "DKS8",
				"Methodics" =>  "ПМ 234-26" ,
				"ProjectsDates" =>  "" ,
				"SampleComment" => "",
				"ProjectsNames" => ""
		)
	 const MEASUREMENTS_SPECIFICATION_GOUP_NAME = Ref("measurements_specification") 
"""
Structure to store the information on temperature experiment

- `sensors_locations` - locations of thermocouples in `m` 
- `material` - material name
- `comments` - comments on sample
- `thickness` - sample total thickness in `m`
"""
	mutable struct SampleProperties
		sensors_locations :: OrderedDict{String , Float64} 
		material::String
		comments::String
		thickness::Float64
		SampleProperties(;
			sensors_locations  = OrderedDict{String , Float64}() , 
			material ::String = "unknown" , 
			comments ::String = "" , 
			thickness ::Float64 = 0.0) = new(sensors_locations , material , comments , thickness)
	end
"""
Structure to connect selected sensors names, sample properties and `WinPosProject`
The idea is that some sensors data should be assiciated with the sample properties in
a single structure. The data (thermocouple measurements) can be trimmed in time and 
combined into single matrix with `[time , T1...TN]` with [`combine_selected_data`](@ref)
Further several [`DataSelector`](@ref) objects can be combined into [`DataSelectorsGroup`](@ref)
- `tmin_tmax` - time range which can be used to select the timerange to trim measured data 
- `selected_names` - stores selected names 
- `sample_properties` - stores sample properties [`SampleProperties`](@ref)

"""
	struct DataSelector{P}

		project :: P
		tmin_tmax :: MVector{2 , Float64}
		selected_names :: OrderedSet{String}
		sample_properties :: SampleProperties

		function DataSelector(;
			project::WinPosProject,
			tmin_tmax::NTuple{2 , Float64} = (-Inf64 , Inf64),
			selected_names::Union{Vector{String} , Nothing} = nothing ,
			sample_properties::Union{SampleProperties , Nothing} = nothing) 

			sample = isnothing(sample_properties) ? SampleProperties() : sample_properties
			_selected_names = isnothing(selected_names) ? OrderedSet{String}() : OrderedSet{String}(filter(Base.Fix1(haskey , project) , selected_names))
			new{typeof(project)}(
				project , 
				MVector(tmin_tmax) , 
				_selected_names , 
				sample
			)
		end
	end
	name(d::DataSelector) = d.project.name
	select!(d::DataSelector , names)  = foreach(names) do n 
		select!(d , n)
	end
	unselect!(d::DataSelector , names ) = foreach(names) do n 
		unselect!(d , n)
	end
	hasname(d::DataSelector , name) = haskey(d.project , name)
	isselected(d::DataSelector , name) = name ∈ d.selected_names
	select!(d::DataSelector , name::String) = hasname(d , name) && push!(d.selected_names , name) 
	select!(d::DataSelector) = foreach(all_names(d)) do name
		push!(d.selected_names , name)
	end
	unselect!(d::DataSelector , name:: String) = isselected(d , name) && delete!(d.selected_names , name)
	unselect!(d::DataSelector) = foreach(WinPos.all_names(d)) do name
		delete!(d.selected_names , name)
	end
	"""
    default_tmin_tmax(d::DataSelector)

Returns the tuple of default values for tmin_tmax range according to the selected data 
"""
	default_tmin_tmax(d::DataSelector) =_default_range(d , Inf , -Inf , :x)

	WinPos.all_names(d::DataSelector) = WinPos.all_names(d.project)

	function fill_data_for_selected!(d::DataSelector)
		for (_ , data_pair_i) in each_selected(d)
			WinPos.fill_data!(data_pair_i)
		end
	end

	"""
    fill_all_data!(d::DataSelector)

Fills all data for associated files 
"""
function fill_all_data!(d::DataSelector)
		WinPos.fill_data!(d.project)
	end

	default_temperature_range(d::DataSelector) = _default_range(d , Inf , -Inf , :y)

	function _default_range(d::DataSelector , tmin , tmax , field_name::Symbol) 
		is_any_selected(d) || return (tmin = tmin , tmax = tmax)
		for (_ , data_pair_i) in each_selected(d)
			WinPos.fill_data!(data_pair_i)
			(_tmin , _tmax)  = extrema(getfield(data_pair_i , field_name))
			(_tmin < tmin) && (tmin = _tmin) 
			(_tmax > tmax) && (tmax = _tmax) 
		end
		return (tmin = tmin , tmax = tmax)
	end
	"""
    selected_names(ds::DataSelector)

Returns selected names in the same order as in the sensors locations ordered dict
"""
function selected_names(ds::DataSelector)
		isempty(ds.selected_names) && return String[] 
		return collect(ds.selected_names)#[k for (k , _) in ds.sample_properties.sensors_locations if isselected(ds , k)]
	end
	is_any_selected(d::DataSelector) = !isempty(d.selected_names)
	tmin(d::DataSelector) = first(d.tmin_tmax)
	tmax(d::DataSelector) = last(d.tmin_tmax)
	thickness(d::DataSelector) = d.sample_properties.thickness
	set_thickness!(d::DataSelector , val::Float64)  = begin 
		@assert val > 0 "Thickness must be greater than zero"
		d.sample_properties.thickness = val
		for (k ,l)  in d.sample_properties.sensors_locations
			if l < 0.0
				d.sample_properties.sensors_locationsp[k] = 0.0
			elseif l > val
				d.sample_properties.sensors_locationsp[k] = val
			end
		end
	end
	set_time_region!(d::DataSelector; tmin = nothing , tmax = nothing) = begin
		isnothing(tmin) || setindex!(d.tmin_tmax ,tmin , 1 )
		isnothing(tmax) || setindex!(d.tmin_tmax ,tmax , 2 )
		return nothing
	end
	"""
    each_selected(d::DataSelector)

Iterator over `DataPairs` which are selected in `DataSelector` returns pair of selected name - `DataPair`
"""
each_selected(d::DataSelector) = Iterators.map(s->Pair(s,d.project[s])  , d.selected_names)  
# Iterators.Filter(sens->isselected(d , last(sens).name)  , d.project) # ( last(sens).name ∈ d.selected_names)


	Base.show(io::IO , p::DataSelector) = begin 
        str = join(string.(selected_names(p)) , " , ")
		locs = join(string.(sensors_locations(p)) , " , ")
        println(io , "  DataSelector over `$(p.project) , selected sensors : $(str) located at $(locs)" )  
    end
	"""
    set_location!(d::DataSelector , name::String , value::Float64)

Function to set the sensor location ( it doesn't sets  name sensor as `selected`)
"""
function set_location!(d::DataSelector , name::String , value::Float64) 
	!hasname(d , name) && error("Inappropriate sensor name ")
	(value < 0.0 || value > thickness(d))  &&  error("Sensor location  must be within the  the sample thickness range 0.0 < $(value) <$(thickness(d)) ")
	push!(d.sample_properties.sensors_locations , name => value)
end
set_location!(d :: DataSelector , name_value_pairs) = foreach(name_value_pairs) do (n , v)
	set_location!(d , n , v)
end
	"""
    sensors_locations(ds :: DataSelector)

returns sensors locations
"""
function sensors_locations(ds :: DataSelector)
		#isempty(ds.sample_properties.sensors_locations) && return Float64[] 
		return map(s->get(ds.sample_properties.sensors_locations,s,-1.0) , selected_names(ds))#[d for (k , d) in ds.sample_properties.sensors_locations if isselected(ds , k)]
	end
		"""
		selected_data_cutted_with_keys(d::DataSelector)

	Returns named tuple with `(data , names)` 
	"""
	selected_data_cutted_with_keys(d::DataSelector) = WinPos.to_matrix(d.project , names = selected_names(d) ,
																tmin = d.tmin_tmax[1] , tmax = d.tmin_tmax[2])

	selected_data_with_keys(d::DataSelector)  = WinPos.to_matrix(d.project , names = selected_names(d)) 


	selected_data(d::DataSelector)  = d |> selected_data_with_keys |> first
	selected_data_cutted(d::DataSelector) = d|>  selected_data_cutted_with_keys |> first

	is_all_locations_assigned(d::DataSelector) = all(Base.Fix1(haskey , d.sample_properties.sensors_locations) , selected_names(d) )
	
	"""
    not_assigned_locations_names(d::DataSelector)

returns the vector of sensors names  which has no assigned locations
"""
function not_assigned_locations_names(d::DataSelector)
		( isempty(d.selected_names) || is_all_locations_assigned(d) ) && return  String[]
		isempty(d.sample_properties.sensors_locations) && return selected_names(d)
		return filter(k -> !haskey(d.sample_properties.sensors_locations , k) , collect( d.selected_names ) )
	end
	"""
    combine_selected_data(d::DataSelector)

returns named tuple with the following data:
 - `time_data` - x-vector for selected signals,
 - `temperature` - matrix of sensors data combined in a single matrix 
 - `initial_distribution` - is just a single value, which can be used as a starting value of sensor's data
 - `sensors_locations` - locations values associated with sensors
"""
function combine_selected_data(d::DataSelector) 
		_sensors_names = selected_names(d)
		(time_data , temperatures)  = WinPos.joindata(d.project , names = _sensors_names ,
										tmin = first(d.tmin_tmax) , tmax = last(d.tmin_tmax) )

		is_all_locations_assigned(d) || error("Before combining assign locations to sensors : $(join(not_assigned_locations_names(d), ',')) ")
		_sensors_locations = sensors_locations(d)

		initial_distribution =  sum(temperatures[1 , :])/size(temperatures , 2)
		return (; 
				time_data = time_data,
				temperatures = temperatures, 
				initial_distribution = initial_distribution,
				sensors_locations = _sensors_locations,
				sample_thickness = d.sample_properties.thickness,
				selected_names = _sensors_names
			)
	end
	struct DataSelectorsGroup
		d::OrderedDict{String , DataSelector}
		name::String
		function DataSelectorsGroup(projects::WinPosProjectsGroup; name::String = "no name")
			d = OrderedDict{String , DataSelector}(
				Pair(n , DataSelector(project = p))  for (n , p) in projects	
			)
			name = (name == "no name") ? projects.name : name
			new(d , name)
		end
		function DataSelectorsGroup(p::Pair{String ,  <:DataSelector} ...)
			n = name(last(first(p)))
			new(OrderedDict{String , DataSelector}(p...) , n)
		end
		function DataSelectorsGroup(p::T ...) where T <: Pair{String , WinPosProject}
			n = WinPos.name(last(first(p)))
			new(OrderedDict{String , DataSelector}( pj[1]=>DataSelector(project = pj[2]) for  pj in p ) , n)
		end
		DataSelectorsGroup(d::OrderedDict{String  , DataSelector}; name = "no name") = new(d , name)

	end
	
	select!(d::DataSelectorsGroup , name::String , selected_names) = haskey(d.d , name) && select!(d.d[name] , selected_names)
	unselect!(d::DataSelectorsGroup , name::String , selected_names) = haskey(d.d , name) && unselect!(d.d[name] , selected_names)
	unselect!(d::DataSelectorsGroup) = foreach(d.d) do (_,di)
		unselect!(di)
	end
	
	selected_data(d :: DataSelectorsGroup) = Iterators.map(selected_data , values(d.d))
	selected_data_cutted(d :: DataSelectorsGroup ) = Iterators.map(selected_data_cutted , values(d.d))
	combine_selected_data(d :: DataSelectorsGroup ) = Iterators.map(combine_selected_data , values(d.d))

	function fun_map(a::Pair , f::Function)
		(k , v)  = (first(a) , last(a))
		return Pair(k , f(v))
	end
	fun_map_wrap(f::Function) = Base.Fix2(fun_map , f)
	# (k , v) -> Pair(k  , selected_data(v))
	selected_data_with_keys(d :: DataSelectorsGroup) = Iterators.map( fun_map_wrap(selected_data_with_keys) , d.d)
	# (k , v) -> Pair(k , selected_data_cutted_with_keys(v))
	selected_data_cutted_with_keys(d :: DataSelectorsGroup ) = Iterators.map( fun_map_wrap(selected_data_cutted_with_keys)  , d.d)
	# (k , v) -> Pair(k ,combine_selected_data(v))
	combine_selected_data_with_keys(d :: DataSelectorsGroup ) = Iterators.map( fun_map_wrap(combine_selected_data)  , d.d)	

	"""
    fill_data_for_selected!(d::DataSelectorsGroup)

Loads data from associated winpos files for all selected sensors for each element of selection group
"""
fill_data_for_selected!(d::DataSelectorsGroup) = foreach(d.d) do (_,di)
		fill_data_for_selected!(di)
	end
"""
    fill_all_data!(d::DataSelectorsGroup)

Fills all data (not only selected), all data! 
"""
fill_all_data!(d::DataSelectorsGroup) = foreach(d.d) do (_,di)
	fill_all_data!(di)
end
	function joined_selected_names(d::DataSelectorsGroup; delim::String = ":")
		j_names = String[]
		for (k_i , d_i) in d.d
			for s_i in d_i.selected_names
				push!(j_names , join((k_i , s_i) , delim ))
			end
		end
		return j_names
	end
	Base.getindex(d::DataSelectorsGroup , key::String) = d.d[key] 
	Base.getindex(d::DataSelectorsGroup , i::Int) = (i <= length(d.d)) ? d.d[iterate(keys(d.d) , i)[1]] : error("out of range")
	Base.keys(d::DataSelectorsGroup) = keys(d.d)

	WinPos.all_names(d::DataSelectorsGroup) = [k => WinPos.all_names(di) for (k , di) in d.d]
	selected_names(d::DataSelectorsGroup) = [k => selected_names(di) for (k , di) in d.d]

	function default_temperature_range(dsg::DataSelectorsGroup)
		(Tmin  , Tmax) = (Inf , -Inf)
		for (k , d_i) in dsg.d
			(_Tmin , _Tmax) = default_temperature_range(d_i)
			(_Tmin < Tmin) && (Tmin = _Tmin)
			(_Tmax > Tmax) && (Tmax = _Tmax)

		end
		return (tmin = Tmin , tmax = Tmax)
	end

	function Base.convert(::Type{WinPos.WinPosProjectsGroup} , d::DataSelectorsGroup)
		d_d = OrderedCollections.OrderedDict{String , WinPos.WinPosProject}()
		paths = Vector{String}()
		for ( _ , d_i ) in d.d
			d_d[d_i.project.name] = d_i.project 
			push!( paths , d_i.project.path )
		end
		compath = common_path(paths)
		return WinPos.WinPosProjectsGroup(basename(compath) , d_d , compath)
	end

	"""
    common_path(paths)

Takes a vector of paths and returns the common part
"""
function common_path(paths)
			split_paths = splitpath.(paths)
			min_len = minimum(length.(split_paths))
			common = String[]
			for i in 1 : min_len
				segment = split_paths[1][i]
				if all(p -> p[i] == segment, split_paths)
					push!(common, segment)
				else
					break
				end
			end
			return joinpath(common...)
	end

	"""
    WinPos.export_to_hdf5(projects::DataSelectorsGroup,
		 	fullfilename::String ; 
        	opentype::String = "w" , 
			overwrite_groups::Bool= true ,
			group_name::Union{String , Nothing}= nothing , 
        	add_path_info::Bool = true)

Exports `DataSelectorsGroup` to HDF5 file
"""
function WinPos.export_to_hdf5(projects::DataSelectorsGroup,
		 	fullfilename::String ; 
        	opentype::String = "w" , 
			overwrite_groups::Bool= true ,
			group_name::Union{String , Nothing}= nothing , 
        	add_path_info::Bool = true)

		root_name = !isnothing(group_name) ? group_name : projects.name 
		#(fullfilename, opentype) = WinPos._check_hdf5_filename_opentype(fullfilename , projects , root_name , opentype)
		h5open(fullfilename, opentype) do h5
			g_root = WinPos._delete_if_overwrite_or_create_group!(h5 , root_name , overwrite_groups)   
			add_path_info && (attributes(g_root)["path"] = fullfilename)
			WinPos._set_if_overwrite_or_create_attribute!(g_root , "type" , "DataSelectorsGroup" , true)

			foreach(projects.d) do (p_name, data_selector)
				proj_branch = WinPos._delete_if_overwrite_or_create_group!(g_root , p_name , overwrite_groups)

				WinPos._set_if_overwrite_or_create_attribute!(proj_branch , "type" , "DataSelector" , true)

				WinPos.add_winpos_proj_to_hdf5!(proj_branch , data_selector.project.name , 
							data_selector.project , 
							overwrite_groups , add_path_info)
				# adding combined_data to the project group
				add_combined_data_to_hdf5!(proj_branch , data_selector)
				add_sample_properties_to_hdf5!(proj_branch , data_selector)
				proj_branch["selected_sensors"] = selected_names(data_selector)
				proj_branch["tmin"] = tmin(data_selector)
				proj_branch["tmax"] = tmax(data_selector)
			end
		end
	end
	function WinPos.load_from_hdf5(file_name::Union{String , HDF5.File , HDF5.Group}  ,::Type{D}) where D <: DataSelectorsGroup
		return load_data_selectors_group_from_hdf5(file_name)
	end
	function load_data_selectors_group_from_hdf5(file :: String)
		@assert isfile(file) "Not a file $(file)"
		return h5open(file) do h5 
			#	root_node_name = first(keys(h5))
    		#	root_node = h5[root_node_name]
			load_data_selectors_group_from_hdf5(h5)
		end
	end
	function load_data_selectors_group_from_hdf5(branch :: Union{HDF5.File , HDF5.Group})
		root_type = WinPos.read_data_type_from_hdf5(branch)
		(isnothing(root_type) || (root_type != "DataSelectorsGroup")) && error("Incorrect branch type")
		root_node_key = first(keys(branch))
		root_node = branch[root_node_key]
		data_selectors_dict = OrderedDict{String , DataSelector}()
		for k in keys(root_node)
			# @show k
			cur_node = root_node[k]
			# @show WinPos.read_data_type_from_hdf5(cur_node)
			WinPos.check_data_type(cur_node , DataSelector) || continue
			proj_name = ""
			for q in keys(cur_node)
				sub_node = cur_node[q]
				WinPos.check_data_type(sub_node , WinPos.WinPosProject) || continue
				 # winpos project in current data selected_sensors
				proj_name = q
				break  
			end
			_wprj = WinPos.load_winpos_project_from_hdf5(cur_node , proj_name)
			_smpl_props = read_sample_properties_from_hdf5(cur_node)
			_selected_names = read(cur_node , "selected_sensors")
			_tmin_tmax = (read(cur_node["tmin"]) , read(cur_node["tmax"]))
			_data_selector = DataSelector(project = _wprj  , selected_names = _selected_names ,
			 	tmin_tmax = _tmin_tmax , sample_properties = _smpl_props)
			data_selectors_dict[k] = _data_selector
		end
		DataSelectorsGroup(data_selectors_dict , name = root_node_key) 
	end
	"""
    add_sample_properties_to_hdf5!(proj_branch , data_selector::DataSelector)

Adds sample properties to `proj_branch` from `DataSelector` object
"""
function add_sample_properties_to_hdf5!(proj_branch , data_selector::DataSelector)
		sample_branch = WinPos._delete_if_overwrite_or_create_group!(proj_branch , "sample_properties" , true)
		for fi in fieldnames(SampleProperties)
			fi_str = String(fi)
			val = getfield(data_selector.sample_properties , fi)
			if  isa(val , AbstractString)
				WinPos._set_if_overwrite_or_create_attribute!(sample_branch , fi_str , val , true)
			elseif isa(val , Number)
				sample_branch[fi_str] = val
			end
			
			sens = WinPos._delete_if_overwrite_or_create_group!(sample_branch , "sensors_locations" , true)
			for (k , v) in data_selector.sample_properties.sensors_locations
				sens[k] = v
			end
		end
	end
	function read_sample_properties_from_hdf5(root_branch )
		sample_properties = SampleProperties()
		smpl_props_branch = haskey(root_branch , "sample_properties") ?  root_branch["sample_properties"] : error("Branch has $(keys(root_branch)) but no `sample_properties`")
		sens_locs_branch = haskey(smpl_props_branch , "sensors_locations") ?  smpl_props_branch["sensors_locations"] : error("Branch has $(keys(root_branch)) but no `sensors_locations`")
		sensors_locations = OrderedDict{String , Float64}()
		for k in keys(sens_locs_branch)
			sensors_locations[k] = read(sens_locs_branch[k])
		end
		sample_properties.sensors_locations = sensors_locations
		for k in keys(smpl_props_branch)
			val  = read(smpl_props_branch[k])
		end
		for f_i in fieldnames(SampleProperties)
			f_i_str = String(f_i)
			haskey(smpl_props_branch , f_i_str)  || continue 
			val = read(smpl_props_branch[f_i_str])
			isa(val , Number) || continue
			setfield!(sample_properties , f_i , val)
		end
		for k in keys(attributes(smpl_props_branch))
			hasfield(SampleProperties , Symbol(k))
			val = read(attributes(smpl_props_branch)[k])
			setfield!(sample_properties , Symbol(k) , val)
		end
		return sample_properties
	end 
	function add_combined_data_to_hdf5!(proj_branch , data_selector::DataSelector)
		combined_data_branch = WinPos._delete_if_overwrite_or_create_group!(proj_branch , "combined_data" , true)
		data_combined = combine_selected_data(data_selector)	
		for (name , val) in pairs(data_combined)
			name_str = String(name)
			if haskey(combined_data_branch , name_str)
				delete_object(combined_data_branch , name_str)
			end
			combined_data_branch[name_str] = val
		end	
	end
	add_measurements_specification_to_hdf5(hdf5::AbstractString ,comments_dict::AbstractDict; specification_node_key = MEASUREMENTS_SPECIFICATION_GOUP_NAME[])= h5open(hdf5, "r+") do file
		add_measurements_specification_to_hdf5(file , comments_dict , specification_node_key = specification_node_key)
	end
	function add_measurements_specification_to_hdf5(hdf5::Union{HDF5.File , HDF5.Group} , comments_dict::AbstractDict; specification_node_key = MEASUREMENTS_SPECIFICATION_GOUP_NAME[])
		 	g = WinPos._delete_if_overwrite_or_create_group!(hdf5 , specification_node_key , true )
		    for (k, v) in comments_dict
				g[k] = v
		    end
	end
	"""
		Initially, all measurements are stored in one folder. However, a single sample may involve
	multiple measurement projects, and a specific material may have several different samples.

	"""
	DataConnector

	function read_measurements_specification(hdf5::Union{HDF5.File , HDF5.Group}; 
							specification_node_key::String = MEASUREMENTS_SPECIFICATION_GOUP_NAME[])

		root = WinPos.find_first_node(hdf5 , specification_node_key )
		isnothing(root) && return nothing 
		ms = OrderedDict{String , Union{String , Vector{String} , Number}}()
		for k in keys(DEFAULT_MEASUREMENTS_SPECIFICATION)
			haskey(root , k) && push!(ms , k => read(root[k]))
		end
		return ms
	end
	function read_measurements_specification(hdf5::AbstractString; 
							specification_node_key::String = MEASUREMENTS_SPECIFICATION_GOUP_NAME[])
		isfile(hdf5) || error("Not  a file $(hdf5)")
		return h5open(hdf5 , "r") do io 
			read_measurements_specification(io , specification_node_key = specification_node_key)
		end
	end
	function read_measurements_specification_string(hdf5_full_file::AbstractString , ; 
							specification_node_key::String = MEASUREMENTS_SPECIFICATION_GOUP_NAME[])
		isfile(hdf5_full_file) || return "not a file"
		key_val_iterator = read_measurements_specification(hdf5_full_file , specification_node_key = specification_node_key)
		return join(["$k : $v" for (k, v) in key_val_iterator], "\n")
	end
	
	function winpos_projects_date_names_to_string(d::DataSelectorsGroup)
		dates = ""
		names = ""
		foreach(d.d) do (_,d) 
			name = WinPos.name(d.project)
			is_first_loop = (length(names) == 0) 
			names *= is_first_loop ? "$(name)" : " , $(name)"
			date_str = !contains(name , "_") ?  "unknown" : first(eachsplit(name , "_")) 
			dates *= is_first_loop ? "$(date_str)" : " , $(date_str)"
			
		end
		return (names , dates)
	end
end
