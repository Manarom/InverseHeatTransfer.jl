#using Pkg
#Pkg.activate("notebook")

module CustomPlutoFunctions
 using PlutoUI , HypertextLiteral


 CustomPlutoFunctions
"""
    multi_values(pluto_ui_slider::Type{T} , 
                    data::AbstractDict, 
                    range_constructor::Union{Function , Nothing} = nothing ; 
                    title = "no title", 
                    common_gui_kwargs...) where T <: Union{RangeSlider , Slider}

Makes bindable gui for pluto from each entry of the all_data dictionary
`range_constructor` is a function which must return (tmin , tmax, steplength) from each 
`data` entry, title is the title of combined ui and `common_gui_kwargs`  name - value pairs 
ransfered directly to `pluto_ui_slider` constructor

"""
function multi_values(pluto_ui_slider::Type{T} , 
                    data::AbstractDict, 
                    range_constructor::Union{Function , Nothing} = nothing ; 
                    title = "no title", 
                    common_gui_kwargs...) where T <: Union{RangeSlider , Slider}

    _opertor = isnothing(range_constructor) ? (_) -> (0.0 , 10.0 , 1000) : range_constructor
	_range_constructor(d) = begin 
		(tmin , tmax , npoints) = _opertor(d)
		range(tmin , tmax , npoints)
	end
	PlutoUI.combine() do Child
		@htl("""
		<h6>$title</h6>
		<ul>
		$([
			@htl("<li>$(k): $(Child( k, pluto_ui_slider(_range_constructor(d); default =_range_constructor(d), common_gui_kwargs...)))</li>")
			for (k , d) in data
		])
		</ul>
		""")
	end
end
PlutoUI.CheckBox(args...; default , kwargs...) = PlutoUI.CheckBox(args...;kwargs...)
"""
    multi_values(pluto_ui_element::Type{T} , names::Union{Tuple , AbstractVector}  ; 
            title = "no title", default_values=nothing , defaults = nothing, common_gui_kwargs... ) where T <: Union{NumberField , TextField , CheckBox}


`default_values` - are values set as input arguments to the ui element
`defaults` - values transfered as default keyword argument
"""
function multi_values(  pluto_ui_element::Type{T} , names::Union{Tuple , AbstractVector}  ; 
                        title = "no title",
                        default_values=nothing , 
                        defaults = nothing, 
                        common_gui_kwargs... 
                    ) where T <: Union{NumberField , TextField , CheckBox}
	    
        if isnothing(default_values) 
            fillv =  pluto_ui_element <: PlutoUI.NumberField ? range(0.0, 10, step=1e-2) : pluto_ui_element <: PlutoUI.CheckBox ? false : 60         
            default_values =  ntuple( length(names)) do _
               return fillv
            end
        end
        if isnothing(defaults)
            fillv2 =  pluto_ui_element <: PlutoUI.NumberField ? 0.0 : pluto_ui_element <: PlutoUI.CheckBox ? false : ""
            defaults = ntuple( length(names) ) do _ 
               return fillv2     
            end
        end
        PlutoUI.combine() do Child
            @htl("""
            <h6>$title</h6>
            <ul>
            $([
                @htl("<li>$(name): $(Child(name, pluto_ui_element(deflt ; default = defdef, common_gui_kwargs...)))</li>")
                for (name , deflt , defdef) in zip( names, default_values , defaults)
            ])
            </ul>
            """)
        end
    end

    """
    paired_selected_names(all_data::AbstractDict , unary_operator::Function)

Function returns tuple of Pairs of dictionary key - result of `unary_operator` call 
on corresponding `value`
"""
function paired_selected_names(all_data::AbstractDict , unary_operator::Function)
        #out = Vector{Pair{String , String}}()
        kv = collect(keys(all_data))
        return ntuple(length(all_data)) do i
            k_i = kv[i]
            v = unary_operator(all_data[k_i])
            return Pair(k_i , v)
        end
    end
    function paired_selected_names(all_data::Union{Tuple , Vector{T}} , unary_operator::Function) where T <: Pair
       kv = collect(keys(all_data))
        return ntuple(length(all_data)) do i
            k_i = kv[i]
            v = unary_operator(last(all_data[k_i]))
            return Pair(k_i , v)
        end
    end
    """
    multi_values(pluto_ui_element::Type{T} , 
                unary_operator::Function ,
                all_data::Union{AbstractDict , Vector{Pair}} ; 
                default_values=nothing , 
                title = "No title" , common_gui_kwargs... ) where T <: Union{NumberField , TextField , CheckBox}


"""
function multi_values(pluto_ui_element::Type{T} , 
                unary_operator::Function ,
                all_data::Union{AbstractDict , Vector{Pair}} ; 
                default_values=nothing , defaults = nothing,
                title = "No title" , common_gui_kwargs... ) where T <: Union{NumberField , TextField , CheckBox}

        if isnothing(default_values) 
            fillv =  pluto_ui_element <: PlutoUI.NumberField ? range(0.0, 10, step=1e-2) : (pluto_ui_element <: PlutoUI.CheckBox ? false : ""  )      
            default_values =  ntuple( length(names)) do _
               return fillv
            end
        end
        if isnothing(defaults)
            fillv2 =  pluto_ui_element <: PlutoUI.NumberField ? 0.0 : pluto_ui_element <: PlutoUI.CheckBox ? false : ""
            defaults = ntuple( length(all_data) ) do _ 
               return fillv2     
            end
        end    
        PlutoUI.combine() do Child
                @htl("""
                <h6>$title</h6>
                <ul>
                $([
                    @htl("<li>
                        
                        $(join( [n , f] , ":")): $(
                            Child( 
                                join( [n , f] , ":") , pluto_ui_element(df ; common_gui_kwargs...)
                                )
                        )
                    </li>")
                    
                    for (df , (n, f) , defdef) in zip(default_values, paired_selected_names(all_data , unary_operator) , defaults)
                ])
                </ul>
                """)
            end
    end

    """
    multi_values(pluto_ui_element::Type{T} , 
        all_data , unary_operator ) where T <: Union{PlutoUI.MultiSelect , PlutoUI.Select}


"""
function multi_values(pluto_ui_element::Type{T} , 
        all_data , unary_operator ) where T <: Union{PlutoUI.MultiSelect , PlutoUI.Select}

        PlutoUI.combine() do Child
            @htl("""
            <ul>
            $([
                @htl("<li>$(name): \n $(
                    Child(
                    name , 
                    pluto_ui_element(data_selector)
                    )
                )
                </li>")
                
                for (name, data_selector) in paired_selected_names(all_data , unary_operator)
            ])
            </ul>
            """)
        end

    end

end