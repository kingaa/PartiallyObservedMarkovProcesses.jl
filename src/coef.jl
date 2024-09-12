export coef, coef!

"""
`coef` extracts the parameters stored in a *PompObject*.
"""
coef = function (
    object::PompObject,
    name::Union{Symbol,Nothing} = nothing,
    names...
        )
    if isnothing(name)
        object.params
    else
        nm = intersect((name,names...),keys(object.params))
        object.params[nm]
    end
end

"""
`coef!` alters, appends, or (optionally) replaces the parameters stored in a *PompObject*.
"""
coef! = function (
    object::PompObject,
    params::Union{NamedTuple,Nothing} = nothing;
    reset::Bool = false
    )
    if reset
        object.params = params
    elseif !isnothing(params)
        existing = keys(object.params)
        new = keys(params)
        old = setdiff(existing,new)
        object.params = (;params[new]...,object.params[old]...)
    end
    nothing
end
