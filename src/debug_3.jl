include(joinpath(@__DIR__,"WinPos.jl"))
data_folder = joinpath(@__DIR__,"..", "test","test_data", "binary_files")
(_, p) = first(WinPos.find_project_pairs(data_folder))
(t,T) = WinPos.joindata(p,names = ("T1","T2","T3","T5"))