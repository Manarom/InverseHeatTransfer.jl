
module DataConnector
	include("WinPos.jl")
	using StaticArrays , OrderedCollections , RecipesBase , Reexport
	using .WinPos
	export SampleProperties , DataSelector , DataSelectorsGroup , combine_selected_data
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

	"""
    set_location!(d::DataSelector , name::String , value::Float64)

Function to set the sensor location ( it doesn't sets  name sensor as `selected`)
"""
set_location!(d::DataSelector , name::String , value::Float64) = hasname(d , name) && push!(d.sample_properties.sensors_locations , name => value)
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
				sensors_locations = _sensors_locations 
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
		function DataSelectorsGroup(p::T ...) where T <: Pair{String , DataSelector}
			n = name(last(first(p)))
			new(OrderedDict{String , DataSelector}(p...) , n)
		end
		function DataSelectorsGroup(p::T ...) where T <: Pair{String , WinPosProject}
			n = WinPos.name(last(first(p)))
			new(OrderedDict{String , DataSelector}( pj[1]=>DataSelector(project = pj[2]) for  pj in p ) , n)
		end
	end
	select!(d::DataSelectorsGroup , name::String , selected_names) = haskey(d.d , name) && select!(d.d[name] , selected_names)
	unselect!(d::DataSelectorsGroup , name::String , selected_names) = haskey(d.d , name) && unselect!(d.d[name] , selected_names)
	selected_data(d :: DataSelectorsGroup) = Iterators.map(selected_data , values(d.d))
	selected_data_cutted(d :: DataSelectorsGroup ) = Iterators.map(selected_data_cutted , values(d.d))
	combine_selected_data(d :: DataSelectorsGroup ) = Iterators.map(combine_selected_data , values(d.d))
	Base.getindex(d::DataSelectorsGroup , key::String) = d.d[key] 
	Base.getindex(d::DataSelectorsGroup , i::Int) = (i <= length(d.d)) ? d.d[iterate(keys(d.d) , i)[1]] : error("out of range")

	function WinPos.export_to_hdf5(projects::DataSelectorsGroup, fullfilename::Union{String , Nothing}=nothing ; 
        opentype::String="w" , overwrite_groups::Bool= true , group_name::Union{String , Nothing}= nothing  , 
        add_path_info::Bool = false)

		root_name = !isnothing(group_name) ? group_name : projects.name 
		(fullfilename, opentype) = WinPos._check_hdf5_filename_opentype(fullfilename , projects , root_name , opentype)
		h5open(fullfilename, opentype) do h5
			g_root = _delete_if_owerwrite_or_create!(h5 , root_name , overwrite_groups)   
			add_path_info && (attributes(g_root)["path"] = projects.path)
			foreach(projects) do (p_name, project)
				add_winpos_proj_to_hdf5!(g_root , p_name , project , overwrite_groups , add_path_info)
			end
		end
	end
	"""
		Initially, all measurements are stored in one folder. However, a single sample may involve
	multiple measurement projects, and a specific material may have several different samples.

	"""
	DataConnector
end
