export times, timezero, obs, states

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
obs(object::AbstractPompObject) = reshape(pomp(object).data,length(object.data),1,1)

"""
    states(object)

`states` extracts the state trajectory of a *PompObject*.
"""
states(object::AbstractPompObject) = begin
    x = pomp(object).states
    if isnothing(x)
        nothing                 # COV_EXCL_LINE
    else
        reshape(x,length(x),1,1)
    end
end

## For internal use (not exported):

obs!(
    object::PompObject,
    data::Array{<:NamedTuple,N},
) where N = begin
    object.data = vec(data)
    nothing                     # COV_EXCL_LINE
end

states!(object::PompObject,states::Array{<:NamedTuple,N}) where N = begin
    object.states = vec(states)
    nothing                     # COV_EXCL_LINE
end

statezero!(object::PompObject,x0::NamedTuple) = begin
    object.x0 = x0
    nothing                     # COV_EXCL_LINE
end

statezero!(object::PompObject,x0::Array{<:NamedTuple,N}) where N = begin
    statezero!(object,x0[1])
end
