
#using .DataConnector
# using Plots

include(joinpath(@__DIR__ , "DataConnector.jl"))
import .DataConnector as DC
using Test
fold = joinpath(@__DIR__ , "binary_files")



ascii_folder = joinpath(@__DIR__, ".."  , ".." , raw"test\test_data\property_inversion_ansys_new") 
#tempdir()
function _test_winpos_groups_equality(w1 , w2 ; common_sensors = ("T1" , "T2" )  , atol = 1e-8)

    for ((_ , p_l) , (_ , p_init)) in zip(w1  , w2)
        for k in  common_sensors
            @test all(isapprox.(p_l[k].x, p_init[k].x , atol=atol))
            @test all(isapprox.(p_l[k].y, p_init[k].y , atol=atol))
        end
    end
    return true
end
# loading data from binary files in winpos type 
wp = DC.WinPos.load_from_winpos_folder(fold)
wp_ascii = DC.WinPos.load_from_ascii_folder(ascii_folder , name_matcher = "Tmeas")

@testset "DataConnector.jl" begin
    temp_folder =  mktempdir(@__DIR__  , cleanup=true)
    for (type,  writer , reader) in zip(("winpos" , "ascii" , "hdf5"),
                                      (DC.WinPos.write_to_winpos_folder , DC.WinPos.write_to_winpos_folder, DC.WinPos.export_to_hdf5 ),
                                      (DC.WinPos.load_from_winpos_folder, DC.WinPos.load_from_winpos_folder , Base.Fix2(DC.WinPos.load_from_hdf5 , DC.WinPosProjectsGroup) )) 
        
        wp2comp = (type == "ascii") ? wp_ascii : wp
        print("Tetsing write/read projects group as $(type) ...")
        group_folder = joinpath(temp_folder , type)
        isdir(group_folder) || mkdir(group_folder)
        (type != "hdf5") || (group_folder = joinpath(temp_folder , type, "group.hdf5"))
        if !(type == "winpos")
            writer(wp2comp , group_folder) # writing data 
        end
        wp_loaded =reader(group_folder) # reading data 
        DC.WinPos.fill_data!(wp_loaded) # reloading data if empty 
        _test_winpos_groups_equality(wp_loaded , wp2comp , common_sensors = ("T1" , "T2" ) , atol=1e-4)
         println("ok")
    end
    data_selector_folder = joinpath(temp_folder , "data_selector_group")

    print("Testing DataSelectorsGroup ...")

    sample_properties = DataConnector.SampleProperties(thickness = 9.3e-3)
    d1 = DC.DataSelector(project = wp_ascii[1] , sample_properties = sample_properties)
    d2 = DC.DataSelector(project = wp_ascii[2], sample_properties = sample_properties)
    d3 = DC.DataSelector(project = wp_ascii[3], sample_properties = sample_properties)

    locs = ntuple(10) do i 
        name = "T$(i)"
        val = 1e-3*(i - 0.7 - 1e-6)
        Pair(name , val)
    end

    for d in (d1,d2,d3)
        DC.select!(d , ("T1" , "T3" , "T5" , "T6" , "T8" , "T10"  ))
        DC.set_location!(d , locs)
    end
    temp_folder = @__DIR__()
    dc_group = DC.DataSelectorsGroup("t1" => d1 , "t2" => d2, "t3" => d3)
    dc_group_file = joinpath(temp_folder, "data_selector_test.hdf5")
    DC.export_to_hdf5(dc_group , dc_group_file , group_name = "test_combined")
    dc_group_loaded = DC.load_from_hdf5(dc_group_file , DC.DataSelectorsGroup)

    for i in eachindex(dc_group)
        _dc = dc_group[i]
        _dc_l = dc_group_loaded[i]
        @test all(DC.selected_names(_dc) .== DC.selected_names(_dc_l))
        for e1 in DC.each_selected(_dc) 
            @test e1[2].name in  _dc.selected_names

        end
    end
    println("ok")
end

dc_group_file = raw"E:\JULIA\JULIA_DEPOT\dev\InverseHeatTransfer.jl\src\data_utils\data_selector_test.hdf5"
dc_group_loaded = DC.load_from_hdf5(dc_group_file , DC.DataSelectorsGroup)

DC.selected_names(dc_group_loaded)
DC.unselect!(dc_group_loaded)
DC.WinPos.all_names(dc_group_loaded)
DC.selected_names(dc_group_loaded)
DC.select!(dc_group_loaded[1] , ("T2" , "T2" , "T3"))
DC.selected_names(dc_group_loaded)
DC.default_tmin_tmax(dc_group_loaded[1])

for (_ , data_pair_i) in DC.each_selected(dc_group_loaded[1])
    @show extrema(data_pair_i.x)

end



    sample_properties = DataConnector.SampleProperties(thickness = 9.3e-3)
    d1 = DC.DataSelector(project = wp_ascii[1] , sample_properties = sample_properties)
    d2 = DC.DataSelector(project = wp_ascii[2], sample_properties = sample_properties)
    d3 = DC.DataSelector(project = wp_ascii[3], sample_properties = sample_properties)



    locs = ntuple(10) do i 
        name = "T$(i)"
        val = 1e-3*(i - 0.7 - 1e-6)
        Pair(name , val)
    end

    for d in (d1,d2,d3)
        DC.select!(d , ("T1" , "T3" , "T5" , "T6" , "T8" , "T10"  ))
        DC.set_location!(d , locs)
    end

    DC.sensors_locations(d1)
    DC.selected_names(d1)
    DC.unselect!(d1 , ("T1" , "T3" , "T5" , "T6" , "T8" , "T10"  ))
    DC.sensors_locations(d1)
    DC.unselect!(d1 , ("T10" , "T5" , "T7"  ))
    DC.set_location!(d1 , "T10" , 0.0001)
    DC.thickness(d1)
    DC.sensors_locations(d1)
    DC.select!(d1 , "T10")



    
    dc_group = DC.DataSelectorsGroup("t1" => d1 , "t2" => d2, "t3" => d3)