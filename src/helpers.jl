"""
    times(object)

`times` extracts the time vector from a *PompObject*.
"""
times(object::AbstractPompObject) = pomp(object).times

"""
    timezero(object)

`timezero` extracts the zero-time (t0) from a *PompObject*.
"""
timezero(object::AbstractPompObject) = pomp(object).t0

"""
    obs(object)

`obs` extracts the vector of observables from a *PompObject*.
"""
obs(object::AbstractPompObject) = pomp(object).obs
obs(object::AbstractArray{<:AbstractPompObject}) = stack(obs.(object))

"""
    init_state(object)

`init_state` extracts the latent state at time t0.
"""
init_state(object::AbstractPompObject) = pomp(object).init_state
init_state(object::AbstractArray{<:AbstractPompObject}) = init_state.(object)

"""
    states(object)

`states` extracts the latent state trajectory of a *PompObject*.
"""
states(object::AbstractPompObject) = pomp(object).states
states(object::AbstractArray{<:AbstractPompObject}) = stack(states.(object))

"""
    coef(object,names...)

`coef` extracts the parameters stored in an *AbstractPompObject* or array thereof.
"""
coef(object::AbstractPompObject) = pomp(object).params

coef(
    object::AbstractPompObject,
    names::Vararg{Symbol},
) = let
    nm = intersect(names,keys(pomp(object).params))
    pomp(object).params[nm]
end

coef(object::AbstractArray{<:AbstractPompObject}) = coef.(object)

coef(
    object::AbstractArray{<:AbstractPompObject},
    names::Vararg{Symbol},
) = map(x -> coef(x,names...), object)

Base.show(io::IO, object::AbstractPompObject) =
    println(io,"<" * string(typeof(object).name.name) * ">")
