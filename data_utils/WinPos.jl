module WinPos
    using Mmap
    using OrderedCollections, Interpolations, Tables
    using DelimitedFiles , RecipesBase
    using HDF5

    export WinPosProject , DataPair , WinPosProjectsGroup 
    export parse_folder_as_winpos_projects , 
            find_winpos_projects , to_matrix , to_table , 
            DataPair , read_winposfile , write_winposfile , export_to_hdf5 , 
            load_winpos_project_from_hdf5 , load_from_hdf5
    
            
    struct DataPair
        name::String
        xfile::String
        yfile::String
        project::String
        x::Vector{Float64}
        y::Vector{Float64}
        DataPair(name , xfile , yfile , project) = new(name , xfile , yfile , project, Float64[], Float64[])
        DataPair(;name , xfile , yfile , project , x::T , y::T) where T <: Vector{Float64} = begin 
            length(x) == length(y) || error("x and y vectors must be of the same size")
            return new(name , xfile , yfile , project, x, y)
        end
    end
    function cut_by_x(dp::DataPair; xmin::Union{Nothing,T} = nothing, xmax::Union{Nothing,T} = nothing) where T<: Number
        xmin = isnothing(xmin) ? -Inf : xmin
        xmax = isnothing(xmax) ? Inf : xmax
        x = dp.x
        w = @. (x >= xmin) & (x <= xmax)
        return (;x = dp.x[w], y = dp.y[w])

    end
    is_data_filled(d::DataPair) = !isempty(d.x) && !isempty(d.y)
# new version is faster 
    function old_read_winposfile!(vec::Vector{T}, filename::String, ::Type{F}  = Float32) where {F <: Number, T <:Number}
            isfile(filename) || return vec
            file_size = filesize(filename) # returns the size of file in bytes
            n = file_size ÷ sizeof(F) # size of the output vector in Float32's
            resize!(vec, n)
            open(filename , "r") do io 
                temp_buffer = Vector{F}(undef, n)
                read!(io, temp_buffer)
                vec .= temp_buffer 
            end
            return vec
        end

    """
    read_winposfile!(vec::Vector{T}, filename::String, ::Type{F}  = Float32) where {F <: Number, T <:Number}

Reads data from winpos format file to a vector of `Float64`, 
`Type{F}` shows the type of data in file 
"""    
    function read_winposfile!(vec::Vector{T}, filename::String, ::Type{F}=Float32) where {F<:Number, T<:Number}
        isfile(filename) || return vec
        open(filename, "r") do io
            n = filesize(io) ÷ sizeof(F)
            data_map = mmap(io, Vector{F} , n) 
            resize!(vec, n)
            vec .= data_map
        end
        return vec
    end
    function write_winposfile(vec::Vector{T}, filename::String, ::Type{F}  = Float32) where {F <: Number, T <:Number}
        open(filename , "w" ) do io 
            for v in vec
                write(io , F(v))
            end
        end
        return nothing 
    end
    read_winposfile(filename) = read_winposfile!(Float64[],filename)
    """
    fill_data!(d::DataPair)

Loads data from attached files 
"""
function fill_data!(d::DataPair)
        read_winposfile!(d.x , d.xfile)
        read_winposfile!(d.y , d.yfile)
        return nothing
    end
    fill_data!(dp::Pair{String , DataPair}) = fill_data!(last(dp))
    is_data_filled(dp::Pair{String , DataPair}) = is_data_filled(last(dp))
    abstract type AbstractWinPosProject{D} end
    struct WinPosProject <: AbstractWinPosProject{DataPair}
        name::String
        data::OrderedDict{String , DataPair}
        path::String
    end
    mutable struct WinPosProjectsGroup <: AbstractWinPosProject{WinPosProject}
        name::String
        data::OrderedDict{String , WinPosProject}
        path::String
    end
    is_data_filled(dp::Pair{String , WinPosProject}) = all(is_data_filled, last(dp))
    is_data_filled(p::AbstractWinPosProject) = all(is_data_filled , p)

    filter_data_names(proj::AbstractWinPosProject, ::Nothing) = keys(proj.data)
    filter_data_names(proj::AbstractWinPosProject, names) = filter(Base.Fix1(haskey,proj.data), names)
    Base.getindex(p::AbstractWinPosProject, k::String)  = p.data[k]
    Base.getindex(p::AbstractWinPosProject , i::Int) = (i <= length(p.data)) ? p.data[iterate(keys(p.data) , i)[1]] : error("out of range")
    name(p::AbstractWinPosProject) = p.name
    Base.iterate(p::AbstractWinPosProject) = iterate(p.data)
    Base.iterate(p::AbstractWinPosProject , j) = iterate(p.data , j)
    
    Base.length(p::AbstractWinPosProject) = Base.length(p.data)
    Base.eltype(::Type{WinPosProjectsGroup}) = Pair{String, WinPosProject}
    Base.eltype(::Type{WinPosProject}) = Pair{String, DataPair}
    Base.keys(p::AbstractWinPosProject) = keys(p.data)
    Base.haskey(p::AbstractWinPosProject , name) = haskey(p.data , name)
    
    Base.show(io::IO , p::WinPosProject) = begin 
        str = join(string.(keys(p.data)) , " , ")
        println(io , " WinPosProject named `$(p.name)` contains : $(str) sensors data" )  
    end
    Base.show(io::IO , p::WinPosProjectsGroup) = begin 
        str = join(string.(keys(p.data)) , " , ")
        println(io , " WinPosProjects group named $(p.name) contains : $(str) WinPos projects" )  
    end
    fill_data!(p::AbstractWinPosProject) = foreach(fill_data! , p)
    fill_data!(p::Pair{String , T}) where T <: AbstractWinPosProject = fill_data!(last(p))

    fill_data!(p::WinPosProject , names) = foreach(p) do (k,d)
        (k ∈ names) && fill_data!(d)
    end
   """
    write_winpos_project(p::WinPosProject , root_folder; new_name::Union{String , Nothing} = nothing)

Writes project creating folder with its name and files `*.x `and `*.dat` corresponding `DataPair`'s
"""
function write_winpos_project(p::WinPosProject , root_folder; new_name::Union{String , Nothing} = nothing)

        !isdir(root_folder) && error("provide folder")
        _name = isnothing(new_name) ? p.name : new_name
        proj_folder = joinpath(root_folder , _name)
        isdir(proj_folder) || mkdir(proj_folder)
        for (k , d) in p
            is_data_filled(d) ||  fill_data!(d)
            d_file_name = joinpath(proj_folder , k*".dat")
            x_file_name = joinpath(proj_folder , k*".x")
            write_winposfile(d.y , d_file_name)
            write_winposfile(d.x , x_file_name)
        end
    end
no_branch_key_error(branch , k) = error("""There is no key `$(k)`, in branch `$(branch)` , available keys: `$(keys(branch))`""")

function load_from_hdf5(h5_file::Union{String , HDF5.File , HDF5.Group}  , D::DataType) #where D <: AbstractWinPosProject
    assert_data_type(h5_file , D)
    return load_winpos_project_from_hdf5(h5_file)
end
assert_data_type(h5_file::Union{String , HDF5.File , HDF5.Group} , ::Type{D}) where D  = check_data_type(h5_file , D) ? nothing : error("Unknown data type $(nameof(D))") 
function check_data_type(h5_file::Union{String , HDF5.File , HDF5.Group} , ::Type{D}) where D
    data_type = read_data_type(h5_file)
    isnothing(data_type)  && return false 
    ( Symbol(data_type) != nameof(D) ) && return false 
    return true
end
function read_data_type(h5_file::AbstractString) 
    !isfile(h5_file) && return nothing
    return h5open(h5_file) do h5
         read_data_type(h5)
    end
end
function read_data_type(h5::Union{HDF5.File , HDF5.Group})
    root_node_name = first(keys(h5))
    root_node = h5[root_node_name]
    return  read_attribute_otherwise(root_node , "type" , nothing)    
end
"""
    load_winpos_project_from_hdf5(h5file_name::String , project_name::Union{String , Nothing} = nothing)

If `project_name` is unspecified - loads `WinPosProject` or `WinPosProjectsGroup` stored in the `h5file_name`
If `project_name` is specified , loads `WinPosProject` specified by name or returns `nothing` if there is not
projects with such name

"""
function load_winpos_project_from_hdf5(h5file_name::String , project_name::Union{String , Nothing} = nothing) 
    
    !isfile(h5file_name) && error("Incorrect filename $(h5file_name)")

    return h5open(h5file_name) do h5 

        root_node_name = first(keys(h5))
        #@show root_node_name
        root_node = h5[root_node_name]
        
        #root_type =  read_attribute_otherwise(root_node , "type" , DomainError(""" 
                                                                                #attribute `type` is not specified """))
        root_type = read_data_type(h5)
        isnothing(root_type) && error("""  attribute `type` is not specified """)
        is_wpp = root_type == "WinPosProject"
        is_wpg = root_type == "WinPosProjectsGroup"
        if !isnothing(project_name) 
            if is_wpp 
                (project_name != root_node_name) && no_branch_key_error(root_node , project_name)
                root_node = h5
            elseif is_wpg && haskey(root_node , project_name)
                is_wpp = true
                is_wpg = false
            else
                no_branch_key_error(root_node , project_name)
            end
            root_node_name = project_name
        else
            is_wpp && (root_node = h5) 
        end  
        return if is_wpg
                    projs = OrderedDict{String , WinPosProject}()
                    for k in keys(root_node)
                        projs[k] = load_winpos_project_from_hdf5(root_node , k)
                    end
                    WinPosProjectsGroup(root_node_name , projs , h5file_name )
                elseif is_wpp
                    load_winpos_project_from_hdf5(root_node , root_node_name)
                else
                    error("Unknown project type ")
                end
    end 
end 

return_call_throw(a) = a
return_call_throw(f::Function) = f()
return_call_throw(a::Exception) = throw(a)


function read_attribute_otherwise(node , att_key , otherwise_return) 
    attr = attributes(node)
    haskey(attr , att_key) && return read(attr , att_key)
    return_call_throw(otherwise_return)
end

function load_winpos_project_from_hdf5(hdf_branch_handle::Union{HDF5.File, HDF5.Group} ,
                                                    name::Union{String , Nothing} = nothing)#)
    
    p_name= if isnothing(name) 
                first(keys(hdf_branch_handle)) 
            else 
                if haskey(hdf_branch_handle , name ) 
                    name 
                else
                    no_branch_key_error(hdf_branch_handle , name)
                end
            end
    root_node = hdf_branch_handle[p_name]

    p_path = read_attribute_otherwise(root_node , "path" , "") # haskey(attributes(root_node), "path") ? read_attribute(root_node ,"path") : ""
    
    project = WinPosProject(p_name , OrderedDict{String , DataPair}() , p_path)
    for s_name  in keys(root_node)
        g_sensor = root_node[s_name]

        x_data = read(g_sensor, "x")
        y_data = read(g_sensor, "y")

        x_file = read_attribute_otherwise(g_sensor , "xfile" , "")  # haskey(attributes(g_sensor), "xfile") ? read(attributes(g_sensor)["xfile"]) : ""
        y_file = read_attribute_otherwise(g_sensor , "yfile" , "")  # haskey(attributes(g_sensor), "yfile") ? read(attributes(g_sensor)["yfile"]) : ""
        
        project.data[s_name] = DataPair(name = s_name , 
                                        x = copy(x_data), y = copy(y_data),
                                        xfile = x_file, yfile = y_file ,
                                        project = p_name )
    end
    return project
end
"""
    add_winpos_proj_to_hdf5!(hdf_branch_handle::Union{HDF5.File, HDF5.Group} , 
                                        project_name , 
                                        project :: WinPosProject , 
                                        overwrite_groups::Bool  , 
                                        add_path_info::Bool=false)

This function adds `WinPosProject` to opened hdf5 file or its branch 
"""
function add_winpos_proj_to_hdf5!(hdf_branch_handle::Union{HDF5.File, HDF5.Group} , 
                                                        project_name , 
                                                        project :: WinPosProject , 
                                                        overwrite_groups::Bool  , 
                                                        add_path_info::Bool=false)
        
        g_proj = _delete_if_overwrite_or_create_group!(hdf_branch_handle, project_name , overwrite_groups)
        add_path_info && _set_if_overwrite_or_create_attribute!(g_proj , "path" , project.path , true)  
        _set_if_overwrite_or_create_attribute!(g_proj , "type" , "WinPosProject" , true)  

        for (s_name, data_pair) in project

            g_sensor = _delete_if_overwrite_or_create_group!(g_proj, s_name , overwrite_groups)

            add_path_info && _set_if_overwrite_or_create_attribute!(g_sensor , "xfile" , data_pair.xfile , true)
            # (attributes(g_sensor)["xfile"] = data_pair.xfile)
            add_path_info && _set_if_overwrite_or_create_attribute!(g_sensor , "yfile" , data_pair.yfile , true)
            # (attributes(g_sensor)["yfile"] = data_pair.yfile)
            
            g_sensor["x"] = data_pair.x
            g_sensor["y"] = data_pair.y
        end

end
"""
    find_first_node(node::Union{HDF5.File, HDF5.Group}, target_name::String)

Finds the first node matching specified `target_name`
"""
function find_first_node(node::Union{HDF5.File, HDF5.Group}, target_name::String)
    if haskey(node, target_name)
        return node[target_name]
    end
    for name in keys(node)
        child = node[name]
        if child isa HDF5.Group
            result = find_node(child, target_name)
            isnothing(result) || return result
        end
    end
    return nothing
end
function add_attributes!(file::String , d::AbstractDict , project_name::Union{String , Nothing} = nothing, overwrite::Bool = true) 
    !isfile(file) && error("Not a file $(file)")
    h5open(file , "r+") do h5 
        node = isnothing(project_name) ? h5 : find_first_node(h5 , project_name)
        isnothing(node) && error(""" There is no node named `$(project_name)`""")
        add_attributes!(node , d , overwrite)
    end
end
function _set_if_overwrite_or_create_attribute!(node , key , val , overwrite)   
    attrs = attributes(node)
    if haskey(attrs, key) 
        overwrite ? delete_attribute(node, key) : error("Attribute is already existent and marked not overwrite")
    end 
    write_attribute(node, key, val)
    return val
end
add_attributes!(branch::Union{HDF5.File, HDF5.Group} , d::AbstractDict , overwrite::Bool = true) = foreach(d) do d_i
    add_attribute!(branch , d_i , overwrite)
end
"""
    add_attribute(branch::Union{HDF5.File, HDF5.Group} , value::Pair{String} , overwrite::Bool = true)

Function to add attribute to branch , if overwrite == true attributes are forced to be overwritten 

"""
function add_attribute!(branch::Union{HDF5.File, HDF5.Group} , value::Pair{String} , overwrite::Bool = true)
    (k , v) = value
    haskey(attributes(branch) , k) && !overwrite && error("Attribute $(k) cannot be overwritten , set overwrite=true")
    _set_if_overwrite_or_create_attribute!(branch , k , v, overwrite)
end

function add_winpos_proj_to_hdf5!(file_name::String, project_name , 
                                    project :: WinPosProject , 
                                    overwrite_groups::Bool  , 
                                    add_path_info::Bool=false)

    @assert isfile(file_name) "There if no $(file_name) file" 
    h5open(file_name) do h5 
        add_winpos_proj_to_hdf5!(h5 , project_name , project , overwrite_groups , add_path_info)
    end
end
"""
    export_to_hdf5(project::WinPosProject, fullfilename::Union{String , Nothing}=nothing ; 
        opentype::String="w" , overwrite_groups::Bool= true , group_name::Union{String , Nothing}= nothing  , 
        add_path_info::Bool = false , autofill::Bool=true)

Saves `WinPosProject` to hdf5 file
Input arguments 

- `project` 
- `fullfilename` 
- `opentype` attribute for [`HDF5.h5open`](@ref)
- `overwrite_groups` if true all variables with the same name will be rewritten , otherwise returns error
- `group_name` custom name for the project if nothing the name of root group will be the same as the name of the project 
- `add_path_info` if true the path will be added to the attributes
- `autofill` if true the project data will be automatically refilled before saving 

"""
function export_to_hdf5(project::WinPosProject, fullfilename::Union{String , Nothing} = nothing ; 
                                            opentype::String = "w" , overwrite_groups::Bool = true , 
                                            group_name::Union{String , Nothing} = nothing  , 
                                            add_path_info::Bool = true , autofill::Bool = true)

    autofill && !is_data_filled(project) && fill_data!(project)
    root_name = !isnothing(group_name) ? group_name : project.name 
    (fullfilename, opentype) = _check_hdf5_filename_opentype(fullfilename , project , root_name , opentype)
    h5open(fullfilename, opentype) do h5
         add_winpos_proj_to_hdf5!(h5 , root_name , project , overwrite_groups , add_path_info)
    end

end
function export_to_hdf5(projects::WinPosProjectsGroup, fullfilename::Union{String , Nothing}=nothing ; 
        opentype::String="w" , overwrite_groups::Bool= true , group_name::Union{String , Nothing}= nothing  , 
        add_path_info::Bool = true , autofill::Bool=true)

    autofill && !is_data_filled(projects) && fill_data!(projects)
    root_name = !isnothing(group_name) ? group_name : projects.name 
    (fullfilename, opentype) = _check_hdf5_filename_opentype(fullfilename , projects , root_name , opentype)
    h5open(fullfilename, opentype) do h5
        g_root = _delete_if_overwrite_or_create_group!(h5 , root_name , overwrite_groups)   
        add_path_info && (attributes(g_root)["path"] = projects.path)
        attributes(g_root)["type"] = "WinPosProjectsGroup"
        foreach(projects) do (p_name , project)
            add_winpos_proj_to_hdf5!(g_root , p_name , project , overwrite_groups , add_path_info)
        end
    end
end
function _check_hdf5_filename_opentype(fullfilename , projects::AbstractWinPosProject , root_name , opentype)
    fullfilename = if isnothing(fullfilename) 
        joinpath(projects.path , root_name*".hdf5") 
    else
        fullfilename
    end
    opentype = if (opentype != "w") &&  !isfile(fullfilename)
         "w"  
    else
        opentype
    end 
    return (fullfilename , opentype)
end

function _delete_if_overwrite_or_create_group!(root , new_name , overwrite)   
# internal function to check if this group is already exist 
    haskey(root , new_name) && overwrite && delete_object(root , new_name)
    return haskey(root , new_name) ? root[new_name] : create_group(root, new_name)
end

    """
    joindata(proj::WinPosProject; names = nothing, tmin = nothing, tmax = nothing)

Function joins data by names in a single matrix , all `y's` are interpolated according the first x in names

```julia
using WinPos
# wp is a WinPosProject obj
(x , y) = joindata(wp , names = ("T1" , "T2") , xmin = 120.0 , xmax = 1200.0 ])
```

"""
function joindata(proj::WinPosProject; names = nothing, tmin = nothing, tmax = nothing)
        names = filter_data_names(proj , names)
        isempty(names) && return ( ; x = Float64[] , y = Float64[], names = String[])
        n1 = first(names)
        (!haskey(proj.data , n1) || isempty(proj.data)) && return Float64[]
        d1 = proj.data[n1]
        isempty(d1.x) && fill_data!(d1)
        x = copy(d1.x)
        if !isnothing(tmin) || !isnothing(tmax)
            !isnothing(tmin) && filter!(t-> t >= tmin, x )
            (!isnothing(tmax) && (tmax > tmin)) && filter!(t -> t <= tmax, x)
        end
        y = Matrix{Float64}(undef , (length(x) , length(names)) )
        for (i , n_i) in enumerate(names)
            haskey(proj.data , n_i) || continue
            d_i = proj.data[n_i]
            !is_data_consistent(d_i) && try 
                                fill_data!(d_i)
                              catch ex
                                @warn ex
                              end
            !is_data_consistent(d_i) && begin
                                    @warn "Data in $(n_i) is inconsisted and ignored"
                                continue
                            end                  
            y[: , i] = isequal(x , d_i.x) ? d_i.y : linear_interpolation(d_i.x , d_i.y , extrapolation_bc=Line()).(x)
        end
        return ( ; x = x , y = y, names = names)    
    end
    is_data_consistent(a,b) = !isempty(a) && !isempty(b) && (length(a) == length(b))
    is_data_consistent(d::DataPair) = is_data_consistent( d.x , d.y )
    """
    to_matrix(proj::WinPosProject; kwargs...)

Extracts matrix from winpos project data  according to the input variables names 
"""
function to_matrix(proj::WinPosProject; kwargs...) 
        out = joindata(proj ; kwargs...)
        mat = hcat(out.x , out.y)
        names = ("x", out.names...)
        return (; data = mat , names = names)
    end
    function to_table(proj::WinPosProject; kwargs...) 
        (data, names) = to_matrix(proj ; kwargs...)
        return Tables.table(data, header = names)
    end
"""
    find_winpos_projects(root_folder::String; include_subfolders::Bool=false)
    
Searches input folder (and if include_subfolders = true  all subfolders) for winpos files
which has `data_name.x` and `data_name.dat` files pair   
"""
    function find_winpos_projects(root_folder::String; include_subfolders::Bool=false)
        @assert isdir(root_folder) "Incorrect folder name $(root_folder)"
        projects = OrderedDict{String, Tuple{OrderedDict{String, DataPair} , String}}()
        for (folder, subdirs, files) in  walkdir(root_folder)
            # Skip if no .x files (not a project)
            x_files = filter(f -> endswith(lowercase(f), ".x"), files)
            isempty(x_files) && continue
            
            # This is a project folder
            project_name = basename(folder)
            project_path = folder
            project_dict = OrderedDict{String, DataPair}()
            
            for x_file in x_files
                
                data_file = replace(x_file, ".x" => ".dat")
                
                # Check if both files exist
                isfile(joinpath(folder, x_file)) || continue
                isfile(joinpath(folder, data_file)) || continue
                
                x_path = joinpath(folder, x_file)
                data_path = joinpath(folder, data_file)
                
                # Key = filename without extension (no folder path)
                data_name = first(splitext(basename(x_file)))
                project_dict[data_name] = DataPair(data_name, x_path , data_path , project_name)
            end
            
            # Only add if we found valid pairs
            !isempty(project_dict) && (projects[project_name] = (project_dict ,  project_path))
            include_subfolders && break
        end
        
        p_dict = OrderedDict{String, WinPosProject}()
        for (n , p) in projects
            p_dict[n] = WinPosProject(n , p[1] , p[2])
        end
        return WinPosProjectsGroup(basename(root_folder) , p_dict , root_folder)
    end
    """
    parse_folder_as_winpos_projects(dir ; name_matcher::String , variable_name::String = "T" , kwargs...)


    Searches the folder dir for subfolders containing files with `name_matcher`
    and interprets the data in these files as `DataPair` taking each column 
    starting from 2 as dependent variables `Y₁...Yₙ` and first column as the 
    independent.

    The structure of `dir` ,ust be the following:
    ```julia
    # for the following structure `\\dir\\proj1\\T_measured.csv`, `\\dir\\proj2\\T_measured.csv`, 
    # `\\dir\\proj3\\T_measured1.csv` each folder contains a single file matching the `name_matcher`
    d = raw"dir"

    WD = parse_folder_as_winpos_projects(d , name_matcher = "T_measured" , variable_name = "T")

    # now WD is an OrderedDict with keys "proj1", "proj2" and "proj3"
    # each element of WD is a `WinPosProject` with `data` - `OrderedDict`
    # with keys "T1".."TN", where each element is the `DataPair` objetc
```

"""
function parse_folder_as_winpos_projects(dir ; name_matcher::String , variable_name::String = "T" , kwargs...)
		projs = OrderedDict{String , WinPosProject}()
		for f in readdir(dir)
			full_f = joinpath(dir , f)
			isdir(full_f) || continue
			(wp , fname) = parse_ascii_files_to_data_pair(full_f , name_matcher , f ,  variable_name , kwargs...)
			!isempty(wp) || continue
			projs[f] = WinPosProject( basename(full_f) , wp , fname )
		end
		return WinPosProjectsGroup(basename(dir) , projs , dir)
    end
    function parse_ascii_files_to_data_pair(fold , 
                name_matcher::String, project_name::String = "" , 
                variable_name::String = "T"; kwargs...)
            @assert isdir(fold) "Must be dir"
            data = Matrix{Float64}(undef , 0,0)
            full_f = ""
            for f in readdir(fold)
                full_f = joinpath(fold , f)
                isfile(full_f) || continue
                contains(f , name_matcher) || continue
                data = readdlm(full_f; kwargs...)
                break
            end
            t = data[: , 1]
            Tnumb = size(data , 2) - 1
            N = size(data , 1)
            data_names = variable_name .* string.(collect(1 : Tnumb))
            d = OrderedDict{String , DataPair}()
            for i in 2 : Tnumb + 1
                name = data_names[ i -  1]
                d[name] = DataPair(name , "" , "" , project_name) # 			DataPair(name , xfile , yfile , project)
                resize!(d[name].x , N)
                copyto!(d[name].x , t)
                resize!(d[name].y , N)
                copyto!(d[name].y , data[: , i])
            end
	    return (d , full_f)
    end


    @recipe function f(d::DataPair)
        label --> d.name
        return (d.x , d.y )
    end
    @recipe function f(wp::WinPosProject; selected_keys = nothing , xmin = nothing, xmax = nothing)
        has_keys = isnothing(selected_keys)
        for (k , d) in wp
            has_keys || k ∈ selected_keys || continue
            (x , y) = cut_by_x(d , xmin = xmin , xmax = xmax)
            @series begin
                label := "$(k)"
                (x , y)
            end
        end    
    end
    """
        WinPos

    A Julia module for managing, processing, and converting measurement data stored 
    in the WinPos binary format (`.x` and `.dat`).

    The module organizes data into a hierarchy where measurement pairs (X and Y coordinates) 
    are grouped into projects, which can be further aggregated into groups based on 
    specific materials or samples.

    # Data Hierarchy
    - `DataPair`: The core unit containing metadata (file paths) and data vectors (`x`, `y`).
    - `WinPosProject`: A collection of named `DataPair` measurements belonging to a single session.
    - `WinPosProjectsGroup`: A container for multiple projects, facilitating data organization 
    where a single **sample** or **material** involves multiple measurement sessions.

    # Main Functions and Types
    - **Data Management**: `DataPair`, `WinPosProject`, `WinPosProjectsGroup`.
    - **I/O Operations**: `read_winposfile`, `write_winposfile`, `export_to_hdf5`.
    - **Parsing**: `parse_folder_as_winpos_projects`, `find_winpos_projects`.
    - **Processing**: `fill_data!`, `cut_by_x`.
    - **Conversion**: `to_matrix`, `to_table`.

    # Key Features
    - **Efficient I/O**: Fast reading/writing of WinPos binary files via buffer reinterpretation.
    - **HDF5 Integration**: Export and import entire project structures while preserving hierarchy and metadata.
    - **Processing Tools**: Lazy loading (`fill_data!`), data slicing (`cut_by_x`), and compatibility with `Tables.jl` and `RecipesBase.jl`.

    # Data Organization Logic
    The module is designed to handle cases where:
    1. All measurements are initially stored in a single folder.
    2. A **particular sample** has several associated **measurement projects**.
    3. A **particular material** consists of several different **samples**.
    """
    WinPos
end


    