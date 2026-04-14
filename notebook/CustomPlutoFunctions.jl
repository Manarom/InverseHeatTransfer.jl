#using Pkg
#Pkg.activate("notebook")

#module CustomPlutoFunctions
#using PlutoUI , HypertextLiteral


    """
    multi_values(names, default_values=nothing)

Creates a set of NumberField elements named according to `names`, default values can be provided  
"""
function multi_values(names, default_values=nothing)
	    isnothing(default_values) && (default_values = zeros(length(names)))
        PlutoUI.combine() do Child
            @htl("""
            <h6>Thicknesses, mm</h6>
            <ul>
            $([
                @htl("<li>$(name): $(Child(name, NumberField(0:1e3:20, default=deflt)))</li>")
                for (name,deflt) in zip( names, default_values)
            ])
            </ul>
            """)
        end
    end


#end