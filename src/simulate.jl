export simulate

"""
    simulate(object; nsim = 1, params, rinit, rprocess, rmeasure, args...)

Simulate the POMP.
Returns an array of *PompObject*s.
At least the `rinit`, `rprocess`, and `rmeasure` basic components, are needed.
"""
simulate(
    object::ValidPompData = nothing;
    nsim::Integer = 1,
    params::Union{P,AbstractVector{P}} = coef(object),
    rinit::Union{Function,Missing} = missing,
    rprocess::Union{PompPlugin,Missing} = missing,
    rmeasure::Union{Function,Missing} = missing,
    args...,
) where {P<:NamedTuple} = let
    try
        params = val_array(params)
        object = pomp(
            object;
            rinit=rinit,
            rprocess=rprocess,
            rmeasure=rmeasure,
            args...,
        )
        [simulate1(object,params[i])
         for i ∈ eachindex(params), j ∈ 1:nsim]
    catch e
        if hasproperty(e,:msg)
            error("in `simulate`: " * e.msg)
        else
            throw(e)            # COV_EXCL_LINE
        end
    end
end

simulate1(
    object::PompObject,
    params::P,
) where {P<:NamedTuple} = let
    params = val_array(params)
    x0 = rinit(object,params=params,nsim=1)
    x = rprocess(object,x0=x0,params=params)
    y = rmeasure(object,x=x,params=params)
    PompObject(
        object.t0,
        object.times,
        object.timevar,
        object.accumvars,
        params[1],x0[1],vec(x),vec(y),
        object.rinit,
        object.rprocess,
        object.rmeasure,
        object.logdmeasure
    )
end

export simulate_array

"""
    simulate_array(object; nsim = 1, params, rinit, rprocess, rmeasure, args...)

Simulate the POMP.
At least the `rinit`, `rprocess`, and `rmeasure` basic components, are needed.
Return an array containing the simulated sample paths.
"""
simulate_array(
    object::ValidPompData = nothing;
    nsim::Integer = 1,
    params::Union{P,AbstractVector{P}} = coef(object),
    rinit::Union{Function,Missing} = missing,
    rprocess::Union{PompPlugin,Missing} = missing,
    rmeasure::Union{Function,Missing} = missing,
    args...,
) where {P<:NamedTuple} = let
    try
        params = val_array(params)
        object = pomp(
            object;
            rinit=rinit,
            rprocess=rprocess,
            rmeasure=rmeasure,
            args...,
        )
        x0 = POMP.rinit(object,params=params,nsim=nsim)
        x = POMP.rprocess(object,x0=x0,params=params)
        y = POMP.rmeasure(object,x=x,params=params)
        timevar = object.timevar
        t = stack(
            fill(
                NamedTuple{(timevar,)}.(POMP.times(object)),
                length(params),
                nsim
            )
        )
        n = length(t)
        map(
            (t,y,x) -> merge(t,y,x),
            t,y,x
        )
    catch e
        if hasproperty(e,:msg)
            error("in `simulate`: " * e.msg)
        else
            throw(e)            # COV_EXCL_LINE
        end
    end
end

simulate(_...) = error("Incorrect call to `simulate`.")
