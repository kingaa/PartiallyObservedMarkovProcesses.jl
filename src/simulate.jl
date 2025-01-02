using Debugger
# using Revise

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
) where {P<:NamedTuple} = let
    @show " TEST 2 "
    # try
        params = val_array(params)
        object = pomp(
            object;
            rinit= rinit, #ismissing(rinit) ? POMP.rinit : rinit,
            rprocess=rprocess,
            rmeasure=rmeasure,
            args...,
        )

        # [simulate1(object,params[i])
        #  for i ∈ eachindex(params), j ∈ 1:nsim]
        #

        ix = Iterators.product( eachindex(params), 1:nsim )


        x0 = (CONFIG.usethreads) ? 
            Threads.map( jk -> POMP.rinit( object, params = params[jk[1]] ), ix  ) :
            map( jk -> POMP.rinit( object, params = params[jk[1]] ), ix )

        # TODO unwrap 
        x = (CONFIG.usethreads) ? 
            Threads.map( jk -> POMP.rprocess( object, x0 = x0[jk[2]], params = params[jk[1]] ), ix ) :
            map( jk -> POMP.rprocess( object, x0 = x0[jk[2]], params = params[jk[1]] ), ix )

        y = (CONFIG.usethreads) ? 
            Threads.map( jk -> POMP.rmeasure( object, x = x[jk...], params = params[jk[1]] ) , ix ) :
            map( jk -> POMP.rmeasure( object, x = x[jk...], params = params[jk[1]] ), ix )

@bp 
       [ 
            PompObject(
                object.t0,
                object.times,
                object.timevar,
                object.accumvars,
                params[i],x0[j][1],vec(x[i,j]),vec(y[i,j]),
                object.rinit,
                object.rprocess,
                object.rmeasure,
                object.logdmeasure
            )
            for (i,j) in ix  
       ] 
    # catch e
    #     if hasproperty(e,:msg)
    #         error("in `simulate`: " * e.msg)
    #     else
    #         throw(e)            # COV_EXCL_LINE
    #     end
    # end
end

# deprecate ?
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
