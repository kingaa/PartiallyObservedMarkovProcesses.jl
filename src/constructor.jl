export pomp

import DataFrames: DataFrame

abstract type AbstractPompObject{T,P,A,X0,X,Y} end

struct PompObject{
    T <: Time,
    P <: Union{<:NamedTuple,Nothing},
    A <: Union{<:NamedTuple,Nothing},
    X0 <: Union{<:NamedTuple,Nothing},
    X <: Union{Vector{<:NamedTuple},Nothing},
    Y <: Union{Vector{<:NamedTuple},Nothing}
    } <: AbstractPompObject{T,P,A,X0,X,Y}
    t0::T
    times::Vector{T}
    accumvars::A
    params::P
    init_state::X0
    states::X
    obs::Y
    rinit::Union{Function,Nothing}
    rprocess::Union{Function,Nothing}
    rmeasure::Union{Function,Nothing}
    logdmeasure::Union{Function,Nothing}
end

## The following type is valid for the `object` in a call to a POMP function.
ValidPompData = Union{
    Nothing,
    Vector{<:NamedTuple},
    DataFrame,
    AbstractPompObject
}

"""
`pomp` is the constructor for the *PompObject* class.

    pomp(
        data;
        t0, times,
        params,
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
- `params`: parameters. A NamedTuple or vector of NamedTuples.
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
    times::Union{T,AbstractVector{T}},
    params::Union{P,Nothing} = nothing,
    accumvars::Union{<:NamedTuple,Nothing} = nothing,
    rinit::Union{Function,Nothing} = nothing,
    rprocess::Union{Function,Nothing} = nothing,
    rmeasure::Union{Function,Nothing} = nothing,
    logdmeasure::Union{Function,Nothing} = nothing,
) where {Y<:NamedTuple,T1<:Time,T<:Time,P<:NamedTuple} = begin
    if T != T1
        error("`t0` and time-vector must have the same elementary type.")
    end
    times = val_array(Vector{T}(times))
    if t0 > times[1]
        error("`t0` must be no later than first observation time.")
    end
    if any(diff(times).<0)
        error("observation times must be nondecreasing.")
    end
    if !isnothing(data) && length(data) != length(times)
        error("data and times must be of the same length.")
    end
    PompObject(
        t0,
        times,
        accumvars,
        params,
        nothing,
        nothing,
        data,
        rinit,
        rprocess,
        rmeasure,
        logdmeasure
    )
end

import DataFrames: DataFrame, select, eachrow
import InvertedIndices: Not

pomp(
    data::DataFrame;
    t0::T,
    times::Symbol,
    args...,
) where {T<:Time} = begin
    time = getproperty(data,times)::Vector{T}
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

"""
    pomp(object::AbstractPompObject; params=missing, accumvars=missing, rinit=missing, rprocess=missing, rmeasure=missing, logdmeasure=missing)

This form returns a modified version of `object`.
Individual basic components can be modified or removed.
The default is to leave them unchanged.
"""
pomp(
    object::AbstractPompObject;
    params::Union{NamedTuple,Nothing,Missing} = missing,
    accumvars::Union{NamedTuple,Nothing,Missing} = missing,
    rinit::Union{Function,Nothing,Missing} = missing,
    rprocess::Union{Function,Nothing,Missing} = missing,
    rmeasure::Union{Function,Nothing,Missing} = missing,
    logdmeasure::Union{Function,Nothing,Missing} = missing,
) = begin
    if ismissing(params)
        params = pomp(object).params
    end
    if ismissing(accumvars)
        accumvars = pomp(object).accumvars
    end
    if ismissing(rinit)
        rinit = pomp(object).rinit
    end
    if ismissing(rprocess)
        rprocess = pomp(object).rprocess
    end
    if ismissing(rmeasure)
        rmeasure = pomp(object).rmeasure
    end
    if ismissing(logdmeasure)
        logdmeasure = pomp(object).logdmeasure
    end
    PompObject(
        pomp(object).t0,
        pomp(object).times,
        accumvars,
        params,
        nothing,
        nothing,
        pomp(object).obs,
        rinit,
        rprocess,
        rmeasure,
        logdmeasure
    )
end
