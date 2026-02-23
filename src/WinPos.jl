module WinPos

    using OrderedCollections, Interpolations, Tables

    struct DataPair
        name::String
        xfile::String
        yfile::String
        project::String
        x::Vector{Float64}
        y::Vector{Float64}
        DataPair(name , xfile , yfile , project) = new(name , xfile , yfile , project, Float64[], Float64[])
    end
    is_data_filled(d::DataPair) = !isempty(d.x) && !isempty(d.y)
    function read_winposfile!(vec::Vector{T}, filename::String, ::Type{F}  = Float32) where {F <: Number, T <:Number}
        bytes = read(filename)
        n = length(bytes) ÷ 4
        resize!(vec, n)
        copyto!(vec, 1, reinterpret(F, bytes), 1, n)
        return vec
    end
    read_winposfile(filename) = read_winposfile!(Float64[],filename)
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
    """
    joindata(proj::WinPosProject)

Function joins data by names in a single 
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
            is_data_consistent(d_i) && try 
                                fill_data!(d_i)
                              catch ex
                                @show ex
                              end
            !is_data_consistent(d_i) && begin
                                    @warn "Data in $(n_i) is inconsisted and ignored"
                                continue
                            end                  
            y[: , i] = isequal(x , d_i.x) ? d_i.y : linear_interpolation(d_i.x , d_i.y , extrapolation_bc=Line()).(x)
        end
        return ( ; x = x , y = y, names = names)    
    end
    is_data_consistent(a,b) = !isempty(a) && !isempty(b) && length(a) == length(b)
    is_data_consistent(d::DataPair) = is_data_consistent( d.x , d.y )
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
    find_project_pairs(root_folder::String;include_subfolders::Bool=false)
    
Searches input folder (and if include_subfolders = true  all subfolders) for winpos files
which has `data_name.x` and `data_name.dat` files pair   
"""
    function find_project_pairs(root_folder::String;include_subfolders::Bool=false)
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
        return p_dict
    end
end


    