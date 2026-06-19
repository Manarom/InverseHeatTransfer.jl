### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 801835a8-aeed-4ee4-bb45-0896eb0b7b5a
begin 
	import Pkg
	Pkg.activate(@__DIR__)
	Pkg.instantiate()
	using  Revise 
	using PlutoUI 
end

# ╔═╡ 93bebb90-4191-11f1-ad50-9f7954455a16
includet(joinpath(@__DIR__ , "CustomPlutoFunctions.jl"))

# ╔═╡ 0ec0805a-5fdc-4da3-b20a-ef2d00dc3d40
PF = Main.CustomPlutoFunctions

# ╔═╡ 127cef9a-5ab5-4c55-a009-691544b525ec
names = ("a" , "b" )

# ╔═╡ af798650-6e90-436a-9119-2d61e2df4fae
PF.multi_values(TextField , names , title = "text field test" , default_values = (30,50) , defaults = ("dfdf" , "sdsfff"))

# ╔═╡ 0f6e4534-3d65-4655-8787-3be670af1101
PF.multi_values(NumberField , names , title = "number field test" , default_values = (1:30 , 2:10) , defaults = (2 , 3))

# ╔═╡ e811f60e-4a11-47f9-9c5e-0b55fb328479
PF.multi_values(PlutoUI.Slider , names ; title = "number field test" , default_values = (1.0:30 , 2.0:10) , defaults = (2 , 3) , show_value = true)

# ╔═╡ f2e04cb3-d16e-4535-b2ba-d1c5462a2b7b
PF.multi_values(PlutoUI.RangeSlider , names ; title = "number field test" , default_values = (1.0:30 , 2.0:10) , defaults = (2:10 , 3:5) , show_value = true)

# ╔═╡ 2f5bf8fa-16a8-47cc-82ac-c3c1a53da126
PF.multi_values(PlutoUI.CheckBox , names ; title = "checkbox field test" , default_values = (false , false) , defaults = (2:10 , 3:5) , show_value = true )

# ╔═╡ 53112351-c4de-4adb-8459-83fe87870910
PF.multi_values_table((MultiSelect , Select , NumberField) , ("a" , "b" ,"c") ; title = "number field test" , default_values = (string.(1:30) , string.(2:10) , 2:10 ) , defaults = (["2"] , ["4"]  , 1) , show_value = true)

# ╔═╡ d11bf65d-1435-4b48-846a-459df5d7f973
PF.multi_values((MultiSelect , Select , Slider , CheckBox)  , ("a" , "b" ,"c" , "d") ; title = "number field test" , default_values = (string.(1:30) , string.(2:10) , 2:10  , false) , defaults = (["2"] , ["4"]  , 1 , false) , show_value = true)

# ╔═╡ a5ab5c42-9ad7-4568-ade9-73e7d17d3586
multitypes  = (TextField, NumberField , MultiSelect , Select , Slider , CheckBox) 

# ╔═╡ 20fa8e6a-e797-40e7-935e-a3938243228b
@bind fff PF.multi_values_table( multitypes, string.(nameof.(multitypes)) ; title = "number field test" , default_values = (50 , 2:10 , string.(1:30) , string.(2:10) , 2:10  , false) , show_value = true )

# ╔═╡ be35660b-16c1-41be-8396-16f780d4054a
fff

# ╔═╡ 686861ce-2dc8-424f-9de5-5f12060a40c7
PF.multi_values_table(MultiSelect , ("a" , "b" ) ; title = "number field test" , default_values = (string.(1:30) , string.(2:10) ) , defaults = (["2"] , ["4"]) , show_value = true , size=5 , fontsize = 18)

# ╔═╡ 7bbade60-8823-46a0-bf0f-d3097c50816c
range_constructor(t::Tuple) = range(t...)

# ╔═╡ 06fc430c-e5b6-4016-aa9e-fc97999ae696
@bind asd PF.multi_values_table(RangeSlider , Dict("a"=>(0,2.3,1000) , "b" =>(0,2.3,1000) )  , range_constructor;  title = "number field test"  )

# ╔═╡ 0167ae0d-19c9-46c6-bd58-b5054d7b2a07
asd

# ╔═╡ Cell order:
# ╠═801835a8-aeed-4ee4-bb45-0896eb0b7b5a
# ╟─93bebb90-4191-11f1-ad50-9f7954455a16
# ╟─0ec0805a-5fdc-4da3-b20a-ef2d00dc3d40
# ╠═127cef9a-5ab5-4c55-a009-691544b525ec
# ╟─af798650-6e90-436a-9119-2d61e2df4fae
# ╟─0f6e4534-3d65-4655-8787-3be670af1101
# ╠═e811f60e-4a11-47f9-9c5e-0b55fb328479
# ╠═f2e04cb3-d16e-4535-b2ba-d1c5462a2b7b
# ╟─2f5bf8fa-16a8-47cc-82ac-c3c1a53da126
# ╟─53112351-c4de-4adb-8459-83fe87870910
# ╠═d11bf65d-1435-4b48-846a-459df5d7f973
# ╠═a5ab5c42-9ad7-4568-ade9-73e7d17d3586
# ╠═20fa8e6a-e797-40e7-935e-a3938243228b
# ╟─be35660b-16c1-41be-8396-16f780d4054a
# ╟─686861ce-2dc8-424f-9de5-5f12060a40c7
# ╠═7bbade60-8823-46a0-bf0f-d3097c50816c
# ╠═06fc430c-e5b6-4016-aa9e-fc97999ae696
# ╠═0167ae0d-19c9-46c6-bd58-b5054d7b2a07
