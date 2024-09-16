export pomp, pomp!

mutable struct PompObject
    data::Union{Vector{NamedTuple},Nothing}
    t0::Real
    time::Vector{Real}
    params::Union{NamedTuple,Nothing}
    rinit::Union{Function,Nothing}
    rmeasure::Union{Function,Nothing}
end

"""
`pomp` is the constructor for the *PompObject* class.
"""
pomp = function (
    data::Union{DataFrame,Nothing};
    t0::Real,
    times::Symbol,
    params::Union{NamedTuple,Nothing} = nothing,
    rinit::Union{Function,Nothing} = nothing,
    rmeasure::Union{Function,Nothing} = nothing
    )
    time = getproperty(data,times)
    if (t0 > time[1])
        error("`t0` cannot be later than `time[1]`.")
    end
    if (any(diff(time).<0))
        error("observation times must be nondecreasing.")
    end
    data = select(data,Not(times))
    PompObject(
        NamedTuple.(eachrow(data)),
        t0,
        time,
        params,
        rinit,
        rmeasure
    )
end

"""
`pomp!` modifies a *PompObject*.
"""
pomp! = function (
    object::PompObject;
    params::Union{NamedTuple,Nothing} = nothing,
    rinit::Union{Function,Nothing} = nothing,
    rmeasure::Union{Function,Nothing} = nothing
    )
    if !isnothing(params) object.params = params end
    if !isnothing(rinit) object.rinit = rinit end
    if !isnothing(rmeasure) object.rmeasure = rmeasure end
    nothing
end
