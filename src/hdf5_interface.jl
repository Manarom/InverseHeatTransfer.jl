const DC = DataConnector
const WP = DataConnector.WinPos

    function WP.export_to_hdf5(inv_problem::SingleInverseProblem,
		 	fullfilename::String ; 
        	opentype::String = "cw" , 
			group_name::String = INVERSE_PROBLEM_HDF5_GROUPNAME[])    

            h5open(fullfilename, opentype) do h5
               group = WP._delete_if_overwrite_or_create_group!(h5 , group_name , true)
			   WP.export_to_hdf5(inv_problem , group) 
            end
    end

  
    function WP.export_to_hdf5(multiproblems::ParallelInverseProblems,
		 	fullfilename::String ; 
        	opentype::String = "cw" , 
			group_name::String = INVERSE_PROBLEM_HDF5_GROUPNAME[] , 
            add_serialized = true)    

            h5open(fullfilename, opentype) do h5
               group = WP._delete_if_overwrite_or_create_group!(h5 , group_name , true)
               for (i,p) in enumerate(multiproblems.problems)
                    _name = "ip_$(i)"
                    group_i = WP._delete_if_overwrite_or_create_group!(group , _name , true)
			        WP.export_to_hdf5(p , group_i , add_serialized = add_serialized) 
               end
            end
    end
    function WP.export_to_hdf5(inv_problem::SingleInverseProblem , 
                                group::HDF5.Group; 
                                add_serialized::Bool = true, 
                                add_all_stats::Bool = true)

        group_setup = WP._delete_if_overwrite_or_create_group!(group , "setup" , true)
        WP.try_write_struct_to_hdf5(group_setup , inv_problem)

        group_setup["alpha"] = inv_problem.α[]
        group_setup["regularization"] = string(inv_problem.regularization)
        group_setup["psi"] = inv_problem.ψ[]
        group_setup["covariance"] = string(inv_problem.covariance)

        group_op = WP._delete_if_overwrite_or_create_group!(group , "optimizable" , true)

        for (k , o) in pairs(inv_problem.optimizable)
            k_str = string(k)
            o_str = string(o)
            group_op_i = WP._delete_if_overwrite_or_create_group!(group_op , k_str , true)
            attrs = attributes(group_op_i)
            attrs["string"] = o_str
            if applicable(coeffs, o)
                group_op_i["coeffs"] = coeffs(o)
            end
        end
        add_direct_problem(inv_problem , group , add_serialized = add_serialized)    
        
       add_serialized && add_serialized_to_hdf5(inv_problem , group , INVERSE_PROBLEM_HDF5_SERIALIZED_GROUPNAME[] )
       add_all_stats && add_all_stats_to_hdf5(inv_problem , group , ALL_STATS_HDF5_GROUPNAME[] )
    end
    function add_serialized_to_hdf5(data , group ::HDF5.Group , name::String = ""  )
            if length(name) == 0
                name = string(uuid4())
            end    
            if haskey(group , name)
                delete_object(group , name)
            end
            io = IOBuffer()
            serialize(io, data) 
            blob = take!(io) 
            group[name] = blob
            return nothing
    end
    add_direct_problem(inv_problem::SingleInverseProblem , group::HDF5.Group; kwargs...) = add_direct_problem(inv_problem.direct_problem , group ;kwargs...)
    function WinPos.export_to_hdf5(d::HeatTransferProblem , fullfilename::AbstractString ;  
        	opentype::String = "cw" , 
			group_name::String = "direct_problem")

            h5open(fullfilename, opentype) do h5
               group = WP._delete_if_overwrite_or_create_group!(h5 , group_name , true)
			   WP.export_to_hdf5(d , group) 
            end
        end
    function add_direct_problem(d::HeatTransferProblem , group::HDF5.Group; 
                        name::String = "direct_problem",  add_serialized::Bool = false )
        
        group_dir = WP._delete_if_overwrite_or_create_group!(group , name , true)
        WP.try_write_struct_to_hdf5(group_dir , d)
        t = collect(trange(d)) 
        x = collect(xrange(d))
        group_dir["time_grid"] = t                       
        group_dir["x_grid"] = x
        group_dir["upper_BC_type"] = string(OneDHeatTransfer.upper_bc_type(d))
        group_dir["upper_BC"] = d.bc_up.(t)

        group_dir["lower_BC_type"] = string(OneDHeatTransfer.lower_bc_type(d))
        group_dir["lower_BC"] = d.bc_dwn.(t)

        _T = Matrix{Float64}(undef , (100 , 2))
        _T1 = @view _T[:,1]
        _T2 = @view _T[:,2]
        _T1 .=  range(extrema(d.T)... , 100)
        @. _T2 = d.L_f(_T1)
        group_dir["thermal_conductivity"] =_T

        @. _T2 = d.C_f(_T1)
        group_dir["heat_capacity"] = _T

        @. _T2 = d.Ld_f(_T1)
        group_dir["thermal_conductivity_derivative"] =_T
        
        add_serialized && add_serialized_to_hdf5(d , group_dir , "direct_problem_serialized" )
    end