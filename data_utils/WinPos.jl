module WinPos

    using OrderedCollections, Interpolations, Tables
    using DelimitedFiles , RecipesBase
    using HDF5

    export WinPosProject , parse_folder_as_winpos_projects , 
            find_project_pairs , to_matrix , to_table , 
            DataPair , read_winposfile , write_winposfile
    
    struct DataPair
        name::String
        xfile::String
        yfile::String
        project::String
        x::Vector{Float64}
        y::Vector{Float64}
        DataPair(name , xfile , yfile , project) = new(name , xfile , yfile , project, Float64[], Float64[])
    end
    function cut_by_x(dp::DataPair; xmin::Union{Nothing,T} = nothing, xmax::Union{Nothing,T} = nothing) where T<: Number
        xmin = isnothing(xmin) ? -Inf : xmin
        xmax = isnothing(xmax) ? Inf : xmax
        x = dp.x
        w = @. (x >= xmin) & (x <= xmax)
        return (;x = dp.x[w], y = dp.y[w])

    end
    is_data_filled(d::DataPair) = !isempty(d.x) && !isempty(d.y)
    """
    read_winposfile!(vec::Vector{T}, filename::String, ::Type{F}  = Float32) where {F <: Number, T <:Number}

Reads data from winpos format file to a vector of `Float64`, 
`Type{F}` shows the type of data in file 
"""
    function read_winposfile!(vec::Vector{T}, filename::String, ::Type{F}  = Float32) where {F <: Number, T <:Number}
            bytes = read(filename) # reads the data to UInt8
            n = length(bytes) ÷ 4
            resize!(vec, n)
            copyto!(vec, 1, reinterpret(F, bytes), 1, n)
            return vec
        end
    function write_winposfile(vec::Vector{T}, filename::String, ::Type{F}  = Float32) where {F <: Number, T <:Number}
        _vec = F.(vec)
        open(filename , "w" ) do io 
            write(io , _vec )
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
    end
    struct WinPosProject
        name::String
        data::OrderedDict{String , DataPair}
    end
    filter_data_names(proj::WinPosProject, ::Nothing) = keys(proj.data)
    filter_data_names(proj::WinPosProject, names) = filter(Base.Fix1(haskey,proj.data), names)
    Base.getindex(p::WinPosProject, k)  = p.data[k]
    Base.iterate(p::WinPosProject) = Base.iterate(p.data)
    Base.iterate(p::WinPosProject,j) = Base.iterate(p.data,j)
    Base.length(p::WinPosProject) = Base.length(p.data)
    Base.show(io::IO , p::WinPosProject) = begin 
        str = join(string.(keys(p.data)) , ",")
        println(io , " WinPosProject named $(p.name) contains : $(str) " )  
    end
   """
    write_winpos_project(p::WinPosProject , root_folder; new_name::Union{String , Nothing} = nothing)

Writes project creating folder with its name and files corresponding `DataPair`'s
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

    mutable struct WinPosProjectsGroup
        name::String
        projects::OrderedDict{String , WinPosProject}
    end
    Base.get_index(wpg::WinPosProjectsGroup , i::Int) =  
    """
    joindata(proj::WinPosProject)

Function joins data by names in a single matrix , all `y's` are interpolated according the first x in names
"""
function joindata(proj::WinPosProject; names = nothing, tmin = nothing, tmax = nothing)
        names = filter_data_names(proj , names)
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
        projects = OrderedDict{String, OrderedDict{String, DataPair}}()
        for (folder, subdirs, files) in  walkdir(root_folder)
            # Skip if no .x files (not a project)
            x_files = filter(f -> endswith(lowercase(f), ".x"), files)
            isempty(x_files) && continue
            
            # This is a project folder
            project_name = basename(folder)
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
            !isempty(project_dict) && (projects[project_name] = project_dict)
            include_subfolders && break
        end
        p_dict = OrderedDict{String, WinPosProject}()
        for (n , p) in projects
            p_dict[n] = WinPosProject(n , p)
        end
        return WinPosProjectsGroup(root_folder , p_dict)
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
			wp = parse_files_to_data_pair(full_f , name_matcher , f ,  variable_name , kwargs...)
			!isempty(wp) || continue
			projs[f] = WinPosProject(f , wp)
		end
		return WinPosProject(dir , projs)
    end
    function parse_files_to_data_pair(fold , 
                name_matcher::String, project_name::String = "" , 
                variable_name::String = "T"; kwargs...)
            @assert isdir(fold) "Must be dir"
            data = []
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
                d[name] = DataPair(name , full_f , full_f , project_name) # 			DataPair(name , xfile , yfile , project)
                resize!(d[name].x , N)
                copyto!(d[name].x , t)
                resize!(d[name].y , N)
                copyto!(d[name].y , data[: , i])
            end
	    return d
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
    Is a package to work with file projects stored according the following scheme:

    `folder_name//A.dat`
    `folder_name//A.x`
    `folder_name//B.dat`
    `folder_name//B.x`

    `folder_name` is the project name , `A.x` and `A.dat`  store `x` and `y` data 
    
    """
    WinPos
end


    