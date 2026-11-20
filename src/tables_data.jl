


function autocorrelation_stats_table(p::AbstractInverseProblem; is_unweighted::Bool = false)
    stats = autocorrelation_analysis(p , is_unweighted=is_unweighted)
    return autocorrelation_stats_table(autocor_stats; is_unweighted=is_unweighted)
end
"""
    autocorrelation_stats_table(autocor_stats; is_unweighted::Bool = false)

Creates args ready for table data creation 
"""
function autocorrelation_stats_table(autocor_stats; is_unweighted::Bool = false)
	N = length(autocor_stats) # number of problems 
	a1 = autocor_stats[1][1]
	fn = filter(f->isa(getfield(a1 , f) , Number) , fieldnames(typeof(a1))) # tests_names 
	M  = length(fn)
	tbl_vect = Any[]
	for (i , a_i) in enumerate(autocor_stats)
		# a_i - problem stats (vector, each element - corresponds to )
 	   	TPnumber = length(a_i) # length of  	
		tbl_i = Matrix{Any}(undef, (TPnumber ,  M + 1) ) # i'th thermocouple table
		for tp_i in 1 : TPnumber #over couples
			a_i_t = a_i[tp_i] #
			row = Any[]
			push!(row, "P$(i) : T$(tp_i)")
			for (j , f_cur) in enumerate(fn)
				  push!(row , getfield(a_i_t, f_cur))
			end
			tbl_i[tp_i , :] = row
		end
		push!(tbl_vect , tbl_i)
	end
	col_names = Vector{String}(undef,  M + 1)
	col_names[1] = "names"
	for i in 1 : M 
		col_names[i + 1] = "$(fn[i]) "
	end	
	tbb = vcat(tbl_vect...)
	return return (data = tbb  , column_labels=col_names ,  title = "Autocor. tests for $(is_unweighted ? "unweighted res." : "weighted res.")") 
end