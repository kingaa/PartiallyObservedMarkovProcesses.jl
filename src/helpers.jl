export times, timezero, obs, init_state, states, coef

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

"""
    init_state(object)

`init_state` extracts the latent state at time t0.
"""
init_state(object::AbstractPompObject) = pomp(object).init_state

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
