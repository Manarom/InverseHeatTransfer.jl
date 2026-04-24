#using Pkg
#Pkg.activate("notebook")

module CustomPlutoFunctions
 using PlutoUI , HypertextLiteral

"""
Several functions to create gui elements for pluto notebooks with multiple element including mixed types
"""
 CustomPlutoFunctions

 const DEFAULT_ARGS_DICT = Dict( 
                    
                    NumberField => range(0.0, 1, 1000),
                    CheckBox => false,
                    TextField => 50,
                    Select => ["a" , "b"],
                    MultiSelect => ["a" , "b"],    
                    Slider => range(0.0, 1, 1000),   
                    RangeSlider => range(0.0, 1, 1000)    

                )
const DEFAULT_DEFAULS_DICT = Dict(

                    NumberField => 1.0,
                    CheckBox => false,
                    TextField => "",
                    Select => "a",
                    MultiSelect => ["a" , "b"],      
                    Slider => 1.0,   
                    RangeSlider => range(0.0, 1, 1000) 

)
HasNoDefaultUI = Union{CheckBox , NumberField}
HasNoShowValueButHasDefault = Union{CheckBox , NumberField , TextField}
"""
    same_kwargs_ui_constructor(T, args...;  kwargs...)

Function intends to ignore kwargs which are not supported by ui
"""
same_kwargs_ui_constructor(T, args...; show_value = nothing,   kwargs...) = T(args...; kwargs...)
same_kwargs_ui_constructor(T, args...; show_value = nothing, default=nothing,  kwargs...) = T(args...; default=default, kwargs...)
same_kwargs_ui_constructor(::Type{T}, args...; show_value = true,  kwargs...) where T <: Slider = T(args...; show_value=show_value, kwargs...)
same_kwargs_ui_constructor(::Type{T}, args...; default=nothing,  kwargs...) where T <: HasNoDefaultUI = T(args...; kwargs...)
same_kwargs_ui_constructor(::Type{T}, args...; default=nothing, show_value = nothing,  kwargs...) where T <: HasNoDefaultUI = T(args...; kwargs...)


#PlutoUI.CheckBox(args...; default, show_value, kwargs...) = PlutoUI.CheckBox(args...;kwargs...)
#PlutoUI.Select(args...; show_value, kwargs...) = PlutoUI.Select(args...;kwargs...)

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
function multi_values(::Type{T} , 
                    data::AbstractDict, 
                    range_constructor::Union{Function , Nothing} = nothing ; 
                    title = "no title", 
                    common_gui_kwargs...) where T <: Union{RangeSlider , Slider}

    _range_constructor = if !isnothing(range_constructor)
        (d) -> range(range_constructor(d)...)
    else
        (_) -> DEFAULT_ARGS_DICT[T]
    end
	PlutoUI.combine() do Child
		@htl("""
		<h6>$title</h6>
		<ul>
		$([
			@htl("<li>$(k): $(Child( k, same_kwargs_ui_constructor(T , _range_constructor(d); default =_range_constructor(d), common_gui_kwargs...)))</li>")
            for (k , d) in data
		])
		</ul>
		""")
	end
end

"""
    multi_values_table(T::Type, names; title, default_values, defaults, common_gui_kwargs...)

Creates a table where the left column contains the parameter name and the right column 
contains the same type of UI element (T) for each name.

`default_values` - values set as positional input arguments to the UI element.
`defaults` - values passed as the `default` keyword argument (initial state).
"""
function multi_values_table(::Type{T}, names::Union{Tuple, AbstractVector}; 
                            title = "No Title", 
                            column_names = ("Name" , "Value"),
                            fontsize::Int=14,
                            default_values = nothing, 
                            defaults = nothing, 
                            common_gui_kwargs...) where T <: Union{NumberField, TextField, CheckBox, Slider, RangeSlider , Select , MultiSelect}
    
    # Setting default input arguments            
    isnothing(default_values) && (default_values = fill(DEFAULT_ARGS_DICT[T], length(names)))
    # Setting default UI state
    isnothing(defaults) && (defaults = fill(DEFAULT_DEFAULS_DICT[T], length(names)))

    PlutoUI.combine() do Child
        @htl("""
        <div class="ui-table-container" style="font-size: $(fontsize)px;">
            <h6 style="margin-bottom: 10px; font-size: 1.1em; color: #222;">$title</h6>
            <table style="width: 100%; border-collapse: collapse; font-family: sans-serif; border: 1px solid #ddd;">
                <thead>
                    <tr style="border-bottom: 2px solid #ddd; background-color: #f5f5f5;">
                        <th style="text-align: left; padding: 10px; width: 35%; color: #333;">$(column_names[1])</th>
                        <th style="text-align: left; padding: 10px; color: #333;">$(column_names[2])</th>
                    </tr>
                </thead>
                <tbody>
                $([
                    @htl("""
                    <tr style="border-bottom: 1px solid #eee;">
                        <td style="padding: 10px; font-weight: 500; color: #444;">$(name)</td>
                        <td style="padding: 10px;">
                            $(Child(name, same_kwargs_ui_constructor(
                                    T, 
                                    deflt; 
                                    default = defdef, 
                                    common_gui_kwargs...
                                )))
                        </td>
                    </tr>
                    """)
                    for (name, deflt, defdef) in zip(names, default_values, defaults)
                ])
                </tbody>
            </table>
        </div>
        """)
    end
end

"""
    multi_values(pluto_ui_element::Type{T} , names::Union{Tuple , AbstractVector}  ; 
            title = "no title", default_values=nothing , defaults = nothing, common_gui_kwargs... ) where T <: Union{NumberField , TextField , CheckBox}


`default_values` - are values set as input arguments to the ui element
`defaults` - values transfered as default keyword argument
"""
function multi_values(  ::Type{T} , names::Union{Tuple , AbstractVector}  ; 
                        title = "no title",
                        default_values=nothing , 
                        defaults = nothing, 
                        common_gui_kwargs... 
                    ) where T <: Union{NumberField , TextField , CheckBox , Slider , RangeSlider , Select , MultiSelect}
        # setting default input arguments            
	    isnothing(default_values) && (default_values = fill(DEFAULT_ARGS_DICT[T] , length(names)))
        # setting default ui state
        isnothing(defaults) && (defaults= fill(DEFAULT_DEFAULS_DICT[T] , length(names)))

        PlutoUI.combine() do Child
            @htl("""
            <h6>$title</h6>
            <ul>
            $([
                @htl("<li>$(name): $(Child(name, same_kwargs_ui_constructor(T , deflt ; default = defdef, common_gui_kwargs...)))</li>")
                for (name , deflt , defdef) in zip( names, default_values , defaults)
            ])
            </ul>
            """)
        end
    end


"""
    multi_values(  pluto_ui_elements::Union{Tuple , AbstractVector} , names::Union{Tuple , AbstractVector}  ; 
                        title = "no title",
                        default_values=nothing , 
                        defaults = nothing, 
                        common_gui_kwargs... 
                    )

Version for the group of multiple joint gui elements
"""
function multi_values(  pluto_ui_elements::Union{Tuple , AbstractVector} , names::Union{Tuple , AbstractVector}  ; 
                        title = "no title",
                        default_values=nothing , 
                        defaults = nothing, 
                        common_gui_kwargs... 
                    ) 

        if length(pluto_ui_elements) == 1
            return multi_values(  pluto_ui_elements[] , names ; title = title ,
                                    default_values = default_values , 
                                    defaults = defaults , common_gui_kwargs...)
        end

	    @assert length(pluto_ui_elements) == length(names) "Number of gui elements must be the same as the number of names"
        
        if isnothing(default_values)    
            default_values =  ntuple( length(names)) do i
               return DEFAULT_ARGS_DICT[pluto_ui_elements[i]]
            end
        else
            @assert length(default_values) == length(names) "Number of `default_values` must be the same as `names` number"
        end

        if isnothing(defaults)
            defaults =  ntuple( length(names)) do i
               return DEFAULT_DEFAULS_DICT[pluto_ui_elements[i]]
            end
        else
            @assert length(defaults) == length(names) "Number of `default_values` must be the same as `names` number"
        end

        PlutoUI.combine() do Child
            @htl("""
            <h6>$title</h6>
            <ul>
            $([
                @htl("<li>$(name): $(Child(name, same_kwargs_ui_constructor(pluto_ui_element , deflt ; default = defdef, common_gui_kwargs...)))</li>")
                for (pluto_ui_element , name , deflt , defdef) in zip(pluto_ui_elements, names, default_values , defaults)
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
function multi_values(::Type{T} , 
                unary_operator::Function ,
                all_data::Union{AbstractDict , Vector{Pair}} ; 
                default_values=nothing , defaults = nothing,
                title = "No title" , common_gui_kwargs... ) where T <: Union{NumberField , TextField , CheckBox }

        # setting default input arguments            
	    isnothing(default_values) && (default_values = fill(DEFAULT_ARGS_DICT[T] , length(all_data)))
        # setting default ui state
        isnothing(defaults) && (defaults= fill(DEFAULT_DEFAULS_DICT[T] , length(all_data)))  
        
        PlutoUI.combine() do Child
                @htl("""
                <h6>$title</h6>
                <ul>
                $([
                    @htl("<li>
                        
                        $(join( [n , f] , ":")): $(
                            Child( 
                                join( [n , f] , ":") , same_kwargs_ui_constructor(T , df ; default = defdef, common_gui_kwargs...)
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
        all_data , unary_operator ) where T <: PlutoUI.MultiSelect


"""
function multi_values(::Type{T} , 
        all_data , unary_operator ) where T <: Union{MultiSelect , Select}

        PlutoUI.combine() do Child
            @htl("""
            <ul>
            $([
                @htl("<li>$(name): \n $(
                    Child(
                    name , 
                        T(data_selector)
                    )
                )
                </li>")
                
                for (name, data_selector) in paired_selected_names(all_data , unary_operator)
            ])
            </ul>
            """)
        end

    end



    """
    multi_values_table(pluto_ui_elements, names; 
                            title = "Settings Table",
                            default_values = nothing, 
                            defaults = nothing, 
                            common_gui_kwargs...)

Combines multiple ui elements into a single table 
"""
function multi_values_table(pluto_ui_elements, names; 
                            title = "Settings Table",
                            column_names = ("Parameter" , "Value"),
                            fontsize::Int = 14,
                            default_values = nothing, 
                            defaults = nothing, 
                            common_gui_kwargs...)
    
    @assert length(pluto_ui_elements) == length(names) "Number of gui elements must be the same as the number of names"
        
    if isnothing(default_values)    
        default_values =  ntuple( length(names)) do i
            return DEFAULT_ARGS_DICT[pluto_ui_elements[i]]
        end
    else
        @assert length(default_values) == length(names) "Number of `default_values` must be the same as `names` number"
    end

    if isnothing(defaults)
        defaults =  ntuple( length(names)) do i
            return DEFAULT_DEFAULS_DICT[pluto_ui_elements[i]]
        end
    else
        @assert length(defaults) == length(names) "Number of `default_values` must be the same as `names` number"
    end

    PlutoUI.combine() do Child
        @htl("""
        <div class="ui-table-container" style="font-size: $(fontsize)px;">
            <h6 style="margin-bottom: 10px;">$title</h6>
            <table style="width: 100%; border-collapse: collapse;">
                <thead>
                    <tr style="border-bottom: 2px solid #ddd;">
                        <th style="text-align: left; padding: 8px;">$(column_names[1])</th>
                        <th style="text-align: left; padding: 8px;"> $(column_names[2])</th>
                    </tr>
                </thead>
                <tbody>
                $([
                    @htl("""
                    <tr style="border-bottom: 1px solid #eee;">
                        <td style="padding: 8px; font-weight: 500;">$(name)</td>
                        <td style="padding: 8px;">
                            $(Child(name, same_kwargs_ui_constructor(
                                    pluto_ui_element, 
                                    deflt; 
                                    default = defdef, 
                                    common_gui_kwargs...
                                )))
                        </td>
                    </tr>
                    """)
                    for (pluto_ui_element, name, deflt, defdef) in zip(pluto_ui_elements, names, default_values, defaults)
                ])
                </tbody>
            </table>
        </div>
        """)
    end
end
end