include("WinPos.jl")
(_, p) = first(WinPos.find_project_pairs(raw"D:\JuliaDepoth\dev\InverseHeatTransfer.jl\test\test_data\binary_files"))
(t,T) = WinPos.joindata(p,("T1","T2","T3","T5"))