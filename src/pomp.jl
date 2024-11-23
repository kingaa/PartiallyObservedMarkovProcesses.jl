export pomp, pomp!

## Time is the type of time
Time = Real

abstract type AbstractPompObject{T<:Time} end

mutable struct PompObject{T} <: AbstractPompObject{T}
    data::Union{Vector{<:NamedTuple},Nothing}
    t0::T
    times::Vector{T}
    accumvars::Union{NamedTuple,Nothing}
    rinit::Union{Function,Nothing}
    rprocess::Union{Function,Nothing}
    rmeasure::Union{Function,Nothing}
    logdmeasure::Union{Function,Nothing}
end

"""
`pomp` is the constructor for the *PompObject* class.

    pomp(
        data;
        t0, times,
        accumvars,
        rinit, rprocess,
        rmeasure, logdmeasure
        )

## Arguments

- `data`: observations.
  The default constructor takes a vector of NamedTuples as data.
  One can also supply a DataFrame.
- `t0`: zero time, t₀.
- `times`: observation times. If `data` is supplied as a DataFrame, `times` should be a Symbol which is the time variable in the DataFrame.
- `accumvars`: a NamedTuple of state variables to be reset (usually to zero) immediately before each simulation stage.
- `rinit`: simulator of the latent-state distribution at t₀.
  This component should be a function that takes parameters and, optionally, `t0`, the initial time.
- `rprocess`: simulator of the latent-state process.
  This component should be a function that takes states, parameters, and current time (`t`) and returns the updated time and state.
- `rmeasure`: simulator of the measurement process.
  This component should be a function that takes states, parameters, and, optionally, `t`, the current time.
- `logdmeasure`: log pdf of the measurement process.
  This component should be a function that takes data, states, parameters, and, optionally, `t`, the current time.
"""
pomp(
    data::Union{Vector{Y},Nothing} = nothing;
    t0::T1,
    times::Union{T,Vector{T}},
    accumvars::Union{NamedTuple,Nothing} = nothing,
    rinit::Union{Function,Nothing} = nothing,
    rprocess::Union{Function,Nothing} = nothing,
    rmeasure::Union{Function,Nothing} = nothing,
    logdmeasure::Union{Function,Nothing} = nothing,
) where {Y<:NamedTuple,T1<:Time,T<:Time} = begin
    if T != T1
        error("`t0` and time-vector must have the same elementary type.")
    end
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
        accumvars,
        rinit,
        rprocess,
        rmeasure,
        logdmeasure
    )
end

"""
Alternatively, one can construct a *PompObject* from a DataFrame.
In this case, `times` should be the Symbol of the time-variable.
"""
pomp(
    data::DataFrame;
    t0::T,
    times::Symbol,
    args...,
) where {T<:Time} = begin
    time = getproperty(data,times)
    data = NamedTuple.(eachrow(select(data,Not(times))))
    pomp(data;t0=t0,times=time,args...)
end

"""
Given an *AbstractPompObject*, `object`,
`pomp(object)` returns the underlying concrete *PompObject*.
Calling `pomp(object, args...)` returns a copy of `object`, modified
according to `args...`.
"""
pomp(object::PompObject) = object

pomp(
    object::AbstractPompObject;
    args...,
) = begin
    obj = deepcopy(object)
    pomp!(obj;args...)
    obj                         # COV_EXCL_LINE
end

"""
`pomp!(object, args...)` modifies the *PompObject* `object` in place.
Individual basic components can be modified, set, or unset.
"""
pomp!(
    object::PompObject;
    accumvars::Union{NamedTuple,Nothing,Missing} = missing,
    rinit::Union{Function,Nothing,Missing} = missing,
    rprocess::Union{Function,Nothing,Missing} = missing,
    rmeasure::Union{Function,Nothing,Missing} = missing,
    logdmeasure::Union{Function,Nothing,Missing} = missing,
) = begin
    if !ismissing(accumvars)
        object.accumvars = accumvars
    end
    if !ismissing(rinit)
        object.rinit = rinit
    end
    if !ismissing(rprocess)
        object.rprocess = rprocess
    end
    if !ismissing(rmeasure)
        object.rmeasure = rmeasure
    end
    if !ismissing(logdmeasure)
        object.logdmeasure = logdmeasure
    end
    object                      # COV_EXCL_LINE
end

pomp!(
    object::AbstractPompObject;
    args...,
) = begin
    pomp!(pomp(object);args...)
    object                      # COV_EXCL_LINE
end
