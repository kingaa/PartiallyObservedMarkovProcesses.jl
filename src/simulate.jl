export simulate

"""
    simulate(object; params, nsim = 1, args...)

`simulate` simulates the POMP.
`args...` can be used to modify or unset fields.
"""
simulate(
    object::AbstractPompObject;
    params::Union{P,AbstractVector{P}} = coef(object),
    nsim::Integer = 1,
    args...,
) where {P<:NamedTuple} = let
    try
        params = val_array(params)
        object = pomp(object;args...)
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

"""
    simulate(; params, nsim = 1, args...)

`args...` can be used to specify the PompObject.
"""
simulate(
    ;params::Union{P,AbstractVector{P}},
    nsim::Integer = 1,
    args...,
) where {P<:NamedTuple} =
    simulate(pomp(;args...),params=params,nsim=nsim)

simulate1(
    object::AbstractPompObject,
    params::P = coef(object),
) where {P<:NamedTuple} = let
    params = val_array(params)
    x0 = rinit(object,params=params,nsim=1)
    x = rprocess(object,x0=x0,params=params)
    y = rmeasure(object,x=x,params=params)
    PompObject(
        object.t0,object.times,
        object.accumvars,
        params[1],x0[1],vec(x),vec(y),
        object.rinit,
        object.rprocess,
        object.rmeasure,
        object.logdmeasure
    )
end
