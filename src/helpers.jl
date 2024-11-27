export times, timezero, obs, states, coef

"""
    times(object)

`times` extracts the time vector of a *PompObject*.
"""
times(object::AbstractPompObject) = pomp(object).times

"""
    timezero(object)

`timezero` extracts the zero-time (t0) of a *PompObject*.
"""
timezero(object::AbstractPompObject) = pomp(object).t0

"""
    obs(object)

`obs` extracts the time vector of a *PompObject*.
"""
obs(object::AbstractPompObject) = pomp(object).obs

"""
    states(object)

`states` extracts the latent state trajectory of a *PompObject*.
"""
states(object::AbstractPompObject) = pomp(object).states

"""
    coef(object)

`coef` extracts the parameter vector of a *PompObject*.
"""
coef(object::AbstractPompObject) = pomp(object).params

Base.show(io::IO, object::AbstractPompObject) =
    println(io,"<" * string(typeof(object)) * ">.")
