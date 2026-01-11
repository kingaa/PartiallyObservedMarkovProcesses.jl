using DataFrames: AbstractDataFrame, Not, select, eachrow

abstract type AbstractPompObject{T,P,A,X0,X,Y,F,U} end
abstract type PompPlugin end

struct PompObject{
    T <: Time,
    P <: NamedTuple,
    A <: Union{<:NamedTuple,Nothing},
    X0 <: Union{<:NamedTuple,Nothing},
    X <: Union{Vector{<:NamedTuple},Nothing},
    Y <: Union{Vector{<:NamedTuple},Nothing},
    F <: Union{PompPlugin,Nothing},
    U <: NamedTuple,
    } <: AbstractPompObject{T,P,A,X0,X,Y,F,U}
    t0::T
    times::Vector{T}
    timevar::Symbol
    accumvars::A
    params::P
    init_state::X0
    states::X
    obs::Y
    rinit::Union{Function,Nothing}
    rprocess::F
    rmeasure::Union{Function,Nothing}
    logdmeasure::Union{Function,Nothing}
    rprior::Union{Function,Nothing}
    logdprior::Union{Function,Nothing}
    userdata::U
    PompObject(
        ;t0,
        times,
        timevar=:time,
        accumvars=(;),
        params=(;),
        init_state=nothing,
        states=nothing,
        obs=nothing,
        rinit=nothing,
        rprocess=nothing,
        rmeasure=nothing,
        logdmeasure=nothing,
        rprior=nothing,
        logdprior=nothing,
        userdata=(;),
    ) = begin
        params = repair(params)
        userdata = repair(userdata)
        new{
            typeof(t0),
            typeof(params),typeof(accumvars),
            typeof(init_state),typeof(states),typeof(obs),
            typeof(rprocess),typeof(userdata),
        }(
            t0,times,timevar,accumvars,
            params,init_state,states,obs,
            rinit,rprocess,rmeasure,logdmeasure,
            rprior,logdprior,
            userdata
        )
    end
    ## The following constructor should only be used internally because
    ## it is not guaranteed to return a valid `PompObject`.
    PompObject(
        object::AbstractPompObject;
        t0 = missing,
        times = missing,
        timevar = missing,
        accumvars = missing,
        params = missing,
        init_state = missing,
        states = missing,
        obs = missing,
        rinit = missing,
        rprocess = missing,
        rmeasure = missing,
        logdmeasure = missing,
        rprior = missing,
        logdprior = missing,
        userdata = missing,
    ) = begin
        if ismissing(t0)
            t0 = pomp(object).t0
        end
        if ismissing(times)
            times = pomp(object).times
        end
        if ismissing(params)
            params = pomp(object).params
        end
        if ismissing(timevar)
            timevar = pomp(object).timevar
        end
        if ismissing(accumvars)
            accumvars = pomp(object).accumvars
        end
        if ismissing(init_state)
            init_state = pomp(object).init_state
        end
        if ismissing(states)
            states = pomp(object).states
        end
        if ismissing(obs)
            obs = pomp(object).obs
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
        if ismissing(rprior)
            rprior = pomp(object).rprior
        end
        if ismissing(logdprior)
            logdprior = pomp(object).logdprior
        end
        if ismissing(userdata)
            userdata = pomp(object).userdata
        end
        params = repair(params)
        userdata = repair(userdata)
        new{
            typeof(t0),
            typeof(params),typeof(accumvars),
            typeof(init_state),typeof(states),typeof(obs),
            typeof(rprocess),typeof(userdata),
        }(
            t0,times,timevar,
            accumvars,params,
            init_state,states,obs,
            rinit,rprocess,
            rmeasure,logdmeasure,
            rprior,logdprior,
            userdata
        )
    end
end

repair(x::NamedTuple) = x
repair(x::Nothing) = (;)

## The following type is valid for the `object` in a call to most package functions.
ValidPompData = Union{
    Nothing,
    Vector{<:NamedTuple},
    AbstractDataFrame,
    AbstractPompObject
}

"""
`pomp` is the constructor for the *PompObject* class.

    pomp(
        data;
        t0, times, timevar,
        params,
        accumvars,
        rinit, rprocess,
        rmeasure, logdmeasure,
        rprior, logdprior,
        userdata
        )

## Arguments

- `data`: observations.
  The default constructor takes a vector of NamedTuples as data.
  One can also supply an AbstractDataFrame.
- `t0`: zero time, t₀.
- `times`: observation times. If `data` is supplied as a DataFrame (or AbstractDataFrame), `times` should be a Symbol which denotes the time variable in the DfataFrame. NB: the types of `times` and `t0` must match, or an error will be generated.
- `timevar`: optional symbol.  Name of the time variable.
- `params`: parameters. A NamedTuple or vector of NamedTuples.
- `accumvars`: a NamedTuple of state variables to be reset (usually to zero) immediately before each simulation stage.
- `rinit`: simulator of the latent-state distribution at t₀.
  This component should be a function that takes parameters and, optionally, `t0`, the initial time.
  It should return a NamedTuple of state variables.
- `rprocess`: simulator of the latent-state process.
  This component should be a plugin (see [`euler`](@ref), [`onestep`](@ref), and [`discrete_time`](@ref)).
- `rmeasure`: simulator of the measurement process.
  This component should be a function that takes states, parameters, and, optionally, `t`, the current time.
  It should return a NamedTuple of observable variables.
- `logdmeasure`: log pdf of the measurement process.
  This component should be a function that takes data, states, parameters, and, optionally, `t`, the current time.
  It should return a scalar.
- `rprior`: simulator of the prior distribution on parameters.
  This component should be a function that takes parameters and returns a NamedTuple of parameters.
  It should return a NamedTuple of parameters.
- `logdprior`: log pdf of the prior distribution on parameters.
  This component should be a function that takes parameters.
  It should return a scalar.
- `userdata`: an optional NamedTuple containing elements that will be furnished to each of the basic model components.
"""
pomp(
    data::Union{Vector{Y},Nothing} = nothing;
    t0::T1,
    times::Union{T,AbstractVector{T}},
    timevar::Symbol = :time,
    params::Union{P,Nothing} = nothing,
    accumvars::Union{<:NamedTuple,Nothing} = nothing,
    rinit::Union{Function,Nothing} = nothing,
    rprocess::Union{<:PompPlugin,Nothing} = nothing,
    rmeasure::Union{Function,Nothing} = nothing,
    logdmeasure::Union{Function,Nothing} = nothing,
    rprior::Union{Function,Nothing} = nothing,
    logdprior::Union{Function,Nothing} = nothing,
    userdata::Union{<:NamedTuple,Nothing} = nothing,
) where {Y<:NamedTuple,T1<:Time,T<:Time,P<:NamedTuple} = begin
    if T != T1
        error("`t0` and time-vector must have the same elementary type.")
    end
    times = val_array(collect(times))
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
        t0=t0,
        times=times,
        timevar=timevar,
        accumvars=accumvars,
        params=params,
        obs=data,
        rinit=rinit,
        rprocess=rprocess,
        rmeasure=rmeasure,
        logdmeasure=logdmeasure,
        rprior=rprior,
        logdprior=logdprior,
        userdata=userdata
    )
end

pomp(
    data::AbstractDataFrame;
    t0::T,
    times::Symbol,
    args...,
) where {T<:Time} = begin
    time = getproperty(data,times)::AbstractVector{T}
    data = NamedTuple.(eachrow(select(data,Not(times))))
    pomp(data;t0=t0,times=time,timevar=times,args...)
end

"""
Given an *AbstractPompObject*, `object`,
`pomp(object)` returns the underlying concrete *PompObject*.
Calling `pomp(object; args...)` returns a copy of `object`, modified
according to the keyword arguments `args...`.
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
    params::Union{<:NamedTuple,Nothing,Missing} = missing,
    timevar::Union{Symbol,Missing} = missing,
    accumvars::Union{<:NamedTuple,Nothing,Missing} = missing,
    rinit::Union{Function,Nothing,Missing} = missing,
    rprocess::Union{PompPlugin,Nothing,Missing} = missing,
    rmeasure::Union{Function,Nothing,Missing} = missing,
    logdmeasure::Union{Function,Nothing,Missing} = missing,
    rprior::Union{Function,Nothing,Missing} = missing,
    logdprior::Union{Function,Nothing,Missing} = missing,
    userdata::Union{<:NamedTuple,Nothing,Missing} = missing,
) = begin
    PompObject(
        object,
        timevar=timevar,
        accumvars=accumvars,
        params=params,
        rinit=rinit,
        rprocess=rprocess,
        rmeasure=rmeasure,
        logdmeasure=logdmeasure,
        rprior=rprior,
        logdprior=logdprior,
        userdata=userdata,
        init_state=nothing,
        states=nothing
    )
end

pomp(_...) = error("Incorrect call to `pomp`.")
