export simulate

"""
    simulate(object; nsim = 1, params, rinit, rprocess, rmeasure, args...)

`simulate` simulates the POMP.  At least the `rinit`, `rprocess`, and `rmeasure` basic components, are needed.
"""
simulate(
    object::ValidPompData = nothing;
    nsim::Integer = 1,
    params::Union{P,AbstractVector{P}} = coef(object),
    rinit::Union{Function,Missing} = missing,
    rprocess::Union{Function,Missing} = missing,
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
        [simulate1(object,params[j])
         for i ∈ 1:nsim, j ∈ eachindex(params)]
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
