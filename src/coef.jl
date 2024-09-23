export coef, coef!

"""
    coef(object,names...)

`coef` extracts the parameters stored in a *PompObject*.
"""
coef(
    object::AbstractPompObject,
    names...,
) = begin
    if length(names)==0
        pomp(object).params
    else
        nm = intersect(names,keys(pomp(object).params))
        pomp(object).params[nm]
    end
end

"""
    coef!(object,params,reset=false)

`coef!` alters, appends, or (optionally) replaces the parameters stored in a *PompObject*.
"""
coef!(
    object::AbstractPompObject,
    params::Union{<:NamedTuple,Nothing} = nothing;
    reset::Bool = false,
) = begin
    if reset
        pomp(object).params = params
    elseif !isnothing(params)
        existing = keys(pomp(object).params)
        new = keys(params)
        old = setdiff(existing,new)
        pomp(object).params = (;params[new]...,object.params[old]...)
    end
    nothing                     # COV_EXCL_LINE
end
