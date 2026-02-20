module WinPos

    using OrderedCollections

    struct DataPair
        name::String
        xfile::String
        yfile::String
        project::String
        x::Vector{Float64}
        y::Vector{Float64}
        DataPair(name , xfile , yfile , project) = new(name , xfile , yfile , project, Vector{Float64}(undef,1), Vector{Float64}(undef,1))
    end

    function read_winposfile!(vec::Vector{T}, filename::String, ::Type{F}  = Float32) where {F <: Number, T <:Number}
        bytes = read(filename)
        n = length(bytes) ÷ 4
        resize!(vec, n)
        copyto!(vec, 1, reinterpret(F, bytes), 1, n)
        return vec
    end
    read_winposfile(filename) = read_winposfile!(Float64[],filename)
    function fill_data!(d::DataPair)
        read_winposfile!(d.x,d.xfile)
        read_winposfile!(d.y,d.yfile)
    end
    struct WinPosProject
        name::String
        data::OrderedDict{String,DataPair}
    end
"""
    find_project_pairs(root_folder::String)

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
    
    return projects
end
end


    