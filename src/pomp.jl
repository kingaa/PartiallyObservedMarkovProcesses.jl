export pomp, pomp!

abstract type AbstractPompObject end

mutable struct PompObject <: AbstractPompObject
    data::Union{Vector{<:NamedTuple},Nothing}
    t0::Real
    times::Vector{Real}
    params::Union{NamedTuple,Nothing}
    rinit::Union{Function,Nothing}
    rmeasure::Union{Function,Nothing}
    rprocess::Union{Function,Nothing}
    x0::Union{NamedTuple,Nothing}
    states::Union{Vector{<:NamedTuple},Nothing}
end

"""
`pomp` is the constructor for the *PompObject* class.
"""
pomp(
    data::DataFrame;
    t0::Real,
    times::Symbol,
    params::Union{NamedTuple,Nothing} = nothing,
    rinit::Union{Function,Nothing} = nothing,
    rmeasure::Union{Function,Nothing} = nothing,
    rprocess::Union{Function,Nothing} = nothing,
) = begin
    time = getproperty(data,times)
    if (t0 > time[1])
        error("`t0` cannot be later than `time[1]`.")
    end
    if (any(diff(time).<0))
        error("observation times must be nondecreasing.")
    end
    data = NamedTuple.(eachrow(select(data,Not(times))))
    PompObject(
        data,
        t0,
        time,
        params,
        rinit,
        rmeasure,
        rprocess,
        nothing,
        nothing
    )
end

"""
`pomp!` modifies a *PompObject* in place.
One can replace or unset individual fields.
"""
pomp!(
    object::PompObject;
    params::Union{NamedTuple,Nothing,Missing} = missing,
    rinit::Union{Function,Nothing,Missing} = missing,
    rmeasure::Union{Function,Nothing,Missing} = missing,
    rprocess::Union{Function,Nothing,Missing} = missing,
) = begin
    if !ismissing(params)
        object.params = params
        object.x0 = nothing
        object.states = nothing
    end
    if !ismissing(rinit)
        object.rinit = rinit
        object.x0 = nothing
        object.states = nothing
    end
    if !ismissing(rmeasure)
        object.rmeasure = rmeasure
        object.x0 = nothing
        object.states = nothing
    end
    if !ismissing(rprocess)
        object.rprocess = rprocess
        object.x0 = nothing
        object.states = nothing
    end
    object                      # COV_EXCL_LINE
end

"""
`pomp` returns a the *PompObject* underlying an *AbstractPompObject*,
potentially with modifications.
If modifications are made, the original is not changed.
"""
pomp(object::PompObject) = object

pomp(
    object::AbstractPompObject;
    args...,
) = pomp!(deepcopy(pomp(object));args...)
