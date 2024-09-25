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

The default constructor takes a vector of NamedTuples as data.
"""
pomp(
    data::Union{Vector{Y},Nothing} = nothing;
    t0::T,
    times::Union{T,Vector{T}},
    rinit::Union{Function,Nothing} = nothing,
    rprocess::Union{Function,Nothing} = nothing,
    rmeasure::Union{Function,Nothing} = nothing,
    dmeasure::Union{Function,Nothing} = nothing,
) where {Y<:NamedTuple,T<:Real} = begin
    times = val_array(times)
    if !isnothing(data) && length(data) != length(times)
        error("data and times must be of the same length.")
    end
    if t0 > times[1]
        error("`t0` must be no later than first observation time.")
    end
    if any(diff(times).<0)
        error("observation times must be nondecreasing.")
    end
    PompObject{T}(
        data,
        t0,
        times,
        rinit,
        rprocess,
        rmeasure,
        dmeasure
    )
end

"""
Alternatively, one can construct a *PompObject* from a DataFrame.
"""
pomp(
    data::DataFrame;
    t0::T,
    times::Symbol,
    args...,
) where {T<:Real} = begin
    time = getproperty(data,times)
    data = NamedTuple.(eachrow(select(data,Not(times))))
    pomp(data;t0=t0,times=time,args...)
end

"""
`pomp` returns a the *PompObject* underlying an *AbstractPompObject*,
potentially with modifications to the basic model components.
If modifications are made, the original is not changed.
"""
pomp(object::PompObject) = object

pomp(
    object::AbstractPompObject;
    args...,
) = pomp!(deepcopy(pomp(object));args...)

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
