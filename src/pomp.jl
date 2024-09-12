export pomp, pomp!

mutable struct PompObject
    data::Union{DataFrame,Nothing}
    t0::Real
    time::Vector{Real}
    params::Union{NamedTuple,Nothing}
    rinit::Union{Function,Nothing}
end

"""
`pomp` is the constructor for the *PompObject* class.
"""
pomp = function (
    data::Union{DataFrame,Nothing};
    t0::Real,
    times::Symbol,
    params::Union{NamedTuple,Nothing} = nothing,
    rinit::Union{Function,Nothing} = default_rinit
    )
    data = sort(data,times)
    time = getproperty(data,times)
    if (t0 > time[1])
        error(""" `t0` cannot be later than `time[1]`""")
    end
    data = DataFrames.select(data,Not(times))
    PompObject(
        data,
        t0,
        time,
        params,
        rinit
    )
end

"""
`pomp!` modifies a *PompObject*.
"""
pomp! = function (
    object::PompObject;
    params::Union{NamedTuple,Nothing} = nothing,
    rinit::Union{Function,Nothing} = nothing
    )
    if !isnothing(params) object.params = params end
    if !isnothing(rinit) object.rinit = rinit end
    nothing
end

default_rinit = function (;_...)
    error(""""rinit" not defined""")
end
