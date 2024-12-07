export euler, discrete_time, onestep

struct EulerPlugin{F<:Function} <: PompPlugin
    stepfun::F
    stepsize::RealTime
end

struct DiscreteTimePlugin{F<:Function,T<:Time} <: PompPlugin
    stepfun::F
    stepsize::T
end

struct OneStepPlugin{F<:Function} <: PompPlugin
    stepfun::F
end

"""
    discrete_time(stepfun; dt = 1)

The function `stepfun` should advance the state by one time unit.
The magnitude of the time unit is `dt`.
"""
discrete_time(
    stepfun::Function;
    dt::Time = 1,
) = DiscreteTimePlugin(stepfun,dt)

discrete_time(_...) = error("Incorrect call to `discrete_time`.")

"""
    euler(stepfun; dt)

The function `stepfun` should advance the state by an arbitrary time-increment.
The time-increment will be chosen so that equal-sized steps of duration at most `dt` are taken over any desired interval.
"""
euler(
    stepfun::Function;
    dt::Time,
) = EulerPlugin(stepfun,RealTime(dt))

euler(_...) = error("Incorrect call to `euler`.")

"""
    onestep(stepfun)

The function `stepfun` will be called once to advance the latent-state process over an interval of arbitrary duration.
"""
onestep(
    stepfun::Function,
) = OneStepPlugin(stepfun)

onestep(_...) = error("Incorrect call to `onestep`.")

rprocess_step(
    p::DiscreteTimePlugin{F,T},
    t0::T,
    tf::T,
    x::X,
    params::NamedTuple,
) where {F<:Function,T<:Time,X<:NamedTuple} = let
    t = t0
    while t < tf
        @inline x = p.stepfun(;t=t,dt=p.stepsize,x...,params...)::X
        t += p.stepsize
    end
    t,x
end

rprocess_step(
    p::EulerPlugin{F},
    t0::T,
    tf::T,
    x::X,
    params::NamedTuple,
) where {F<:Function,T<:Time,X<:NamedTuple} = let
    n = ceil(Int64,(tf-t0)/p.stepsize)
    if n > 0
        tstep = (tf-t0)/n
        t = t0
        for _ in 1:n
            @inline x = p.stepfun(;t=t,dt=tstep,x...,params...)::X
            t += tstep
        end
    end
    tf,x
end

rprocess_step(
    p::OneStepPlugin{F},
    t0::T,
    tf::T,
    x::X,
    params::NamedTuple,
) where {F<:Function,T<:Time,X<:NamedTuple} = let
    t = t0
    @inline x = p.stepfun(;t=t,dt=tf-t0,x...,params...)::X
    tf,x
end
