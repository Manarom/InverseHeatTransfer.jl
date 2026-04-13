
module DataConnector
	include("WinPos.jl")
	using StaticArrays , OrderedCollections , RecipesBase , Reexport
	using HDF5
	@reexport using .WinPos
	export SampleProperties , DataSelector , DataSelectorsGroup ,
	 combine_selected_data 
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
	unselect!(d::DataSelector , name:: String) = isselected(d , name) && delete!(d.selected_names , name)
	selected_names(d::DataSelector) = collect(d.selected_names)
	tmin(d::DataSelector) = first(d.tmin_tmax)
	tmax(d::DataSelector) = last(d.tmin_tmax)
	thickness(d::DataSelector) = d.sample_properties.thickness
	set_time_region(d::DataSelector; tmin = nothing , tmax = nothing) = begin
		isnothing(tmin) || setindex!(d.tmin_tmax ,tmin , 1 )
		isnothing(tmax) || setindex!(d.tmin_tmax ,tmax , 2 )
		return nothing
	end
	each_selected(d::DataSelector) = Iterators.Filter(sens->( last(sens).name ∈ d.selected_names) , d.project)

	Base.show(io::IO , p::DataSelector) = begin 
        str = join(string.(collect(p.selected_names)) , " , ")
		locs = join(string.(sensors_locations(p)) , " , ")
        println(io , "  DataSelector over `$(p.project) , selected sensors : $(str) located at $(locs)" )  
    end
	"""
    set_location!(d::DataSelector , name::String , value::Float64)

Function to set the sensor location ( it doesn't sets  name sensor as `selected`)
"""
function set_location!(d::DataSelector , name::String , value::Float64) 
	!hasname(d , name) && error("Inappropriate sensor name ")
	(value < 0.0 || value > thickness(d))  &&  error("Sensor location  must be within the  the sample 0.0 < $(value) <$(thickness(d)) ")
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
		isempty(ds.sample_properties.sensors_locations) && return Float64[] 
		return [d for (k , d) in ds.sample_properties.sensors_locations if k ∈ ds.selected_names]
	end
	selected_data(d::DataSelector)  = WinPos.to_matrix(d.project , names = d.selected_names) |> first
	selected_data_cutted(d::DataSelector) = WinPos.to_matrix(d.project , names = d.selected_names ,
															tmin = d.tmin_tmax[1] , tmax = d.tmin_tmax[2]) |> first
	is_all_locations_assigned(d::DataSelector) = all(Base.Fix1(haskey , d.sample_properties.sensors_locations) , d.selected_names )
	
	"""
    not_assigned_locations_names(d::DataSelector)

returns the vector of sensors names  which has no assigned locations
"""
function not_assigned_locations_names(d::DataSelector)
		( isempty(d.selected_names) || is_all_locations_assigned(d) ) && return  String[]
		isempty(d.sample_properties.sensors_locations) && return collect(d.selected_names)
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

		(time_data , temperatures)  = WinPos.joindata(d.project , names = d.selected_names ,
										tmin = first(d.tmin_tmax) , tmax = last(d.tmin_tmax) )

		is_all_locations_assigned(d) || error("Before combining assign locations to sensors : $(join(not_assigned_locations_names(d), ',')) ")
		_sensors_locations = sensors_locations(d)
		initial_distribution =  sum(temperatures[1 , :])/size(temperatures , 2)
		return (; 
				time_data = time_data,
				temperatures = temperatures, 
				initial_distribution = initial_distribution,
				sensors_locations = _sensors_locations,
				sample_thickness = d.sample_properties.thickness
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
	selected_data(d :: DataSelectorsGroup) = Iterators.map(selected_data , values(d.d))
	selected_data_cutted(d :: DataSelectorsGroup ) = Iterators.map(selected_data_cutted , values(d.d))
	combine_selected_data(d :: DataSelectorsGroup ) = Iterators.map(combine_selected_data , values(d.d))
	Base.getindex(d::DataSelectorsGroup , key::String) = d.d[key] 
	Base.getindex(d::DataSelectorsGroup , i::Int) = (i <= length(d.d)) ? d.d[iterate(keys(d.d) , i)[1]] : error("out of range")
	Base.keys(d::DataSelectorsGroup) = keys(d.d)

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
    			#root_node = h5[root_node_name]
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
	"""
		Initially, all measurements are stored in one folder. However, a single sample may involve
	multiple measurement projects, and a specific material may have several different samples.

	"""
	DataConnector
end
