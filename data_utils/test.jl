
#using .DataConnector
using Plots

include(joinpath(@__DIR__() , "WinPos.jl"))

fold = joinpath(@__DIR__ , "binary_files")

ascii_folder = joinpath(@__DIR__ , ".." , raw"test\test_data\property_inversion_ansys_new") 

wp = WinPos.find_winpos_projects(fold)
wp_ascii = WinPos.parse_folder_as_winpos_projects(ascii_folder , name_matcher = "Tmeas")
#d = DataConnector.DataSelector(project = wp[1])



#d_group = DataConnector.DataSelectorsGroup(collect(p for p in wp)...)

WinPos.export_to_hdf5(wp[1])
WinPos.export_to_hdf5(wp)
WinPos.export_to_hdf5(wp_ascii)

proj_file = joinpath(fold , "04.12.2025_RBSN_9d3_full_T1400_10Ks" , "04.12.2025_RBSN_9d3_full_T1400_10Ks.hdf5")
proj_group_file =  joinpath(fold , "binary_files.hdf5" )

wpss = WinPos.load_winpos_project_from_hdf5(proj_file)
wpsss = WinPos.load_winpos_project_from_hdf5(proj_group_file)

wp2 = WinPos.load_winpos_project_from_hdf5(proj_group_file , "04.12.2025_RBSN_9d3_full_T1400_10Ks")

wp3 = WinPos.load_winpos_project_from_hdf5(proj_file , "04.12.2025_RBSN_9d3_full_T1400_10Ks")

plot(wp2)
plot!(wp3)

sum(wp3["T1"].y .- wp2["T1"].y)

WinPos.add_attributes!(proj_file , Dict("asdf" => rand(10)) , "04.12.2025_RBSN_9d3_full_T1400_10Ks")

dd_dd  = WinPos.load_from_hdf5(proj_group_file , WinPos.WinPosProjectsGroup)
##
include("DataConnector.jl")

wp1_data = DataConnector.WinPos.load_winpos_project_from_hdf5(proj_group_file)



d1 = DataConnector.DataSelector(project = wp1_data[1])
d2 = DataConnector.DataSelector(project = wp1_data[2])
d3 = DataConnector.DataSelector(project = wp1_data[3])
for d in (d1,d2,d3)
    DataConnector.select!(d , ("T1" , "T2" ))
    DataConnector.set_location!(d , ("T2" => 10.01, "T3" => 6.7 , "T1" => 3.4))
end
(x , y , names) = DataConnector.combine_selected_data(d1)

ddd = DataConnector.DataSelectorsGroup("t1" => d1 , "t2" => d2, "t3" => d3)
DataConnector.export_to_hdf5(ddd , joinpath(@__DIR__ , "data_selector.hdf5") , group_name = "test_combined")

## 
wp_ascii = DataConnector.WinPos.parse_folder_as_winpos_projects(ascii_folder , name_matcher = "Tmeas")
dd_2 = DataConnector.DataSelectorsGroup(wp_ascii)


for (k , d) in dd_2.d
    DataConnector.select!(d , ("T1" , "T4"  , "T7"))
    DataConnector.set_location!(d , ("T2" => 10.01, "T3" => 6.7 , "T1" => 3.4))
end
##
include("DataConnector.jl")
DataConnector.WinPos.read_data_type(joinpath(@__DIR__ , "data_selector.hdf5"))

dd = DataConnector.WinPos.load_from_hdf5(joinpath(@__DIR__ , "data_selector.hdf5"), DataConnector.DataSelectorsGroup)