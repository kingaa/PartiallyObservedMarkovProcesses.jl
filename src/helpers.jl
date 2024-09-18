export times, timezero, obs

"""
    times(object)

`times` extracts the time vector of a *PompObject*.
"""
times(object::PompObject) = object.times
times(object::AbstractPompObject) = times(pomp(object))

"""
    timezero(object)

`timezero` extracts the zero-time (t0) of a *PompObject*.
"""
timezero(object::PompObject) = object.t0
timezero(object::AbstractPompObject) = timezero(pomp(object))

"""
    obs(object)

`obs` extracts the time vector of a *PompObject*.
"""
obs(object::PompObject) = object.data
obs(object::AbstractPompObject) = obs(pomp(object))

obs!(
    object::AbstractPompObject,
    data::Vector{<:NamedTuple},
) = begin
    pompobj = pomp(object)
    pompobj.data = data
    nothing                     # COV_EXCL_LINE
end
