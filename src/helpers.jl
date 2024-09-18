export times, timezero, obs

"""
    times(object)

`times` extracts the time vector of a *PompObject*.
"""
times(object::PompObject) = object.times

"""
    timezero(object)

`timezero` extracts the zero-time (t0) of a *PompObject*.
"""
timezero(object::PompObject) = object.t0

"""
    obs(object)

`obs` extracts the time vector of a *PompObject*.
"""
obs(object::PompObject) = object.data
