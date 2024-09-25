export times, timezero, obs

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
obs(object::AbstractPompObject) = pomp(object).data

Base.show(io::IO,object::AbstractPompObject) =
    println(io,"object of type <" * string(typeof(object)) * ">.")
