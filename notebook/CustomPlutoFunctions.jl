#using Pkg
#Pkg.activate("notebook")

module CustomPlutoFunctions
 using PlutoUI , HypertextLiteral
    export multi_values , multi_values_table
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

const SupportedUI = Union{NumberField, TextField, CheckBox, Slider, RangeSlider, Select, MultiSelect}

unsupported_kwarg_names(T) = (:show_value,)
unsupported_kwarg_names(::Type{<:Union{Slider, RangeSlider, Scrubbable}}) = (:default,)
unsupported_kwarg_names(::Type{T}) where T <:Union{CheckBox } = (:show_value, :default)

@inline function same_kwargs_ui_constructor(::Type{T}, args...; kwargs...) where T <: SupportedUI
    clean_kwargs = Base.structdiff(values(kwargs), NamedTuple{unsupported_kwarg_names(T)})
    return T(args...; clean_kwargs...)
end


function _prepare_defaults(pluto_ui_elements , names , default_values , defaults)

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

        return (default_values , defaults)
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
function multi_values(  pluto_ui_elements::Union{Tuple , AbstractVector} , 
                        names::Union{Tuple , AbstractVector}  ; 
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

        (_default_values , _defaults) = _prepare_defaults(pluto_ui_elements , names , default_values , defaults)

        PlutoUI.combine() do Child
            @htl("""
            <h6>$title</h6>
            <ul>
            $([
                @htl("<li>$(name): $(Child(name, same_kwargs_ui_constructor(pluto_ui_element , deflt ; default = defdef, common_gui_kwargs...)))</li>")
                for (pluto_ui_element , name , deflt , defdef) in zip(pluto_ui_elements, names, _default_values , _defaults)
            ])
            </ul>
            """)
        end
    end



    """
        multi_values(pluto_ui_element::Type{T} , names::Union{Tuple , AbstractVector}  ; 
                title = "no title", default_values=nothing , defaults = nothing, common_gui_kwargs... ) where T <: Union{NumberField , TextField , CheckBox}


    `default_values` - are values set as input arguments to the ui element
    `defaults` - values transfered as default keyword argument
    """
    function multi_values(  T::Type{<: SupportedUI} , names::Union{Tuple , AbstractVector}  ; 
                            title = "no title",
                            default_values=nothing , 
                            defaults = nothing, 
                            common_gui_kwargs... 
                        )
                uis= ntuple(Returns(T), length(names))
                multi_values(uis , names ; title = title , default_values = default_values , defaults = defaults, common_gui_kwargs... )
        end

        function multi_values(pluto_elements::Union{ Type{<: SupportedUI} , Tuple , AbstractVector} ,
                                names::AbstractDict, 
                                unary_operator = t->t ;
                                title = "no title",
                                defaults = nothing, 
                                common_gui_kwargs... 
                            )
                _names = (keys(names)...,)
                vs = map(unary_operator, (values(names)...,))
                multi_values(pluto_elements , _names ; default_values = vs, defaults = defaults, title = title , common_gui_kwargs...)
        end

   
 
        """
        multi_values_table(pluto_ui_elements, names; 
                                title = "Settings Table",
                                default_values = nothing, 
                                defaults = nothing, 
                                common_gui_kwargs...)

    Combines multiple ui elements into a single table 
    """
    function multi_values_table(pluto_ui_elements::Union{Tuple , AbstractVector} , 
                                names::Union{Tuple , AbstractVector}  ;  
                                title = "Settings Table",
                                column_names = ("Parameter" , "Value"),
                                fontsize::Int = 14,
                                default_values = nothing, 
                                defaults = nothing, 
                                common_gui_kwargs...)
        
        (_default_values , _defaults) = _prepare_defaults(pluto_ui_elements , names , default_values , defaults)

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
                        for (pluto_ui_element, name, deflt, defdef) in zip(pluto_ui_elements, names, _default_values, _defaults)
                    ])
                    </tbody>
                </table>
            </div>
            """)
        end
    end
    function multi_values_table(::Type{T}, names::Union{Tuple, AbstractVector}; 
                                title = "No Title", 
                                column_names = ("Name" , "Value"),
                                fontsize::Int=14,
                                default_values = nothing, 
                                defaults = nothing, 
                                common_gui_kwargs...) where T <: SupportedUI

        uis= ntuple(Returns(T), length(names))
        multi_values_table(uis , names ; title = title , 
                    column_names = column_names , fontsize=fontsize ,
                    default_values = default_values , defaults = defaults ,
                    common_gui_kwargs...)
    end
    function multi_values_table(pluto_elements::Union{ Type{<: SupportedUI} , Tuple , AbstractVector} ,
                                names::AbstractDict, 
                                unary_operator = t->t ;
                                title = "no title",
                                column_names = ("Name" , "Value"),
                                fontsize::Int=14,
                                defaults = nothing, 
                                common_gui_kwargs... 
                            )

                _names = (keys(names)...,)
                vs = map(unary_operator, (values(names)...,))
                multi_values_table(pluto_elements , _names ; column_names = column_names , fontsize = fontsize , default_values = vs, defaults = defaults, title = title , common_gui_kwargs...)
        end

end