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
    rinit::Union{Function,Nothing,Missing} = missing,
    rprocess::Union{PompPlugin,Nothing,Missing} = missing,
    rmeasure::Union{Function,Nothing,Missing} = missing,
    args...,
) where {P<:NamedTuple} = begin
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
end

simulate1(
    object::PompObject,
    params::P,
) where {P<:NamedTuple} = begin
    params = val_array(params)
    x0 = rinit(object,params=params,nsim=1)
    x = rprocess(object,x0=x0,params=params)
    y = rmeasure(object,x=x,params=params)
    PompObject(
        object,
        params=params[1],
        init_state=x0[1],
        states=vec(x),
        obs=vec(y)
    )
end


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
    rinit::Union{Function,Nothing,Missing} = missing,
    rprocess::Union{PompPlugin,Nothing,Missing} = missing,
    rmeasure::Union{Function,Nothing,Missing} = missing,
    args...,
) where {P<:NamedTuple} = begin
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
    map(merge,t,y,x)
end

simulate(_...) = error("Incorrect call to `simulate`.")
