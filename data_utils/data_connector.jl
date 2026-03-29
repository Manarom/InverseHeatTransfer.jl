include("WinPos.jl")
mutable struct DataConnector{P}
	all_names::Vector{String}
	selected_names::Vector{String}
	data :: Matrix{Float64}
	data_cutted :: Matrix{Float64}
	locations :: Vector{Float64}
	total_thickness :: Float64
	project:: P
end
collect_names(p::WinPos.WinPosProject) = collect(keys(p.data))
