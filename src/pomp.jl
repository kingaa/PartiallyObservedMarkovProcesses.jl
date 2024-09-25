export pomp, pomp!

## T is the type of time (T <: Real)
abstract type AbstractPompObject{T} end

mutable struct PompObject{T} <: AbstractPompObject{T}
    data::Union{Vector{<:NamedTuple},Nothing}
    t0::T
    times::Vector{T}
    rinit::Union{Function,Nothing}
    rprocess::Union{Function,Nothing}
    rmeasure::Union{Function,Nothing}
    dmeasure::Union{Function,Nothing}
end

"""
`pomp` is the constructor for the *PompObject* class.
"""
pomp(
    data::DataFrame;
    t0::T,
    times::Symbol,
    rinit::Union{Function,Nothing} = nothing,
    rprocess::Union{Function,Nothing} = nothing,
    rmeasure::Union{Function,Nothing} = nothing,
    dmeasure::Union{Function,Nothing} = nothing,
) where {T<:Real} = begin
    time = getproperty(data,times)
    if (t0 > time[1])
        error("`t0` must be no later than first observation time.")
    end
    if (any(diff(time).<0))
        error("observation times must be nondecreasing.")
    end
    data = NamedTuple.(eachrow(select(data,Not(times))))
    PompObject{T}(
        data,
        t0,
        time,
        rinit,
        rprocess,
        rmeasure,
        dmeasure
    )
end

pomp(
    ;t0::T,
    times::Union{T,Vector{T}},
    rinit::Union{Function,Nothing} = nothing,
    rprocess::Union{Function,Nothing} = nothing,
    rmeasure::Union{Function,Nothing} = nothing,
    dmeasure::Union{Function,Nothing} = nothing,
) where {T<:Real} = begin
    times = val_array(times)
    if (t0 > times[1])
        error("`t0` must be no later than first observation time.")
    end
    if (any(diff(times).<0))
        error("observation times must be nondecreasing.")
    end
    PompObject{T}(
        nothing,
        t0,
        times,
        rinit,
        rprocess,
        rmeasure,
        dmeasure
    )
end

"""
`pomp!` modifies a *PompObject* in place.
One can replace or unset individual fields.
"""
pomp!(
    object::PompObject;
    rinit::Union{Function,Nothing,Missing} = missing,
    rprocess::Union{Function,Nothing,Missing} = missing,
    rmeasure::Union{Function,Nothing,Missing} = missing,
    dmeasure::Union{Function,Nothing,Missing} = missing,
) = begin
    if !ismissing(rinit)
        object.rinit = rinit
    end
    if !ismissing(rprocess)
        object.rprocess = rprocess
    end
    if !ismissing(rmeasure)
        object.rmeasure = rmeasure
    end
    if !ismissing(dmeasure)
        object.dmeasure = dmeasure
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
