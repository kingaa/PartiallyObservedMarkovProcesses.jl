"""
    traj_match_objfun(object, estimvars; rinit, rprocess,
        logdmeasure, logdprior, args...)

Returns a function that compares a model trajectory with the data,
returning minus the sum of the log likelihood and log prior density.
Trajectory matching consists of the minimization of such an objective function.
At least the `rinit`, `rprocess`, and `logdmeasure` basic components are
needed. `args...` can be used to modify or unset additional fields.
`rprocess` must be a `VectorfieldPlugin` (i.e., constructed via a call to [`vectorfield`](@ref)).

The returned function takes as input a vector or tuple of length equal to that
of `estimvars`. This vector is associated, element-for-element, with the symbols in `estimvars`. By default, `estimvars = keys(coef(object))`.
"""
traj_match_objfun(
    object::ValidPompData,
    estimvars::Union{NTuple{N,Symbol},Vector{Symbol},Missing} = missing;
    rinit::Union{Function,Nothing,Missing} = missing,
    rprocess::Union{VectorfieldPlugin,Nothing,Missing} = missing,
    logdmeasure::Union{Function,Nothing,Missing} = missing,
    logdprior::Union{Function,Nothing,Missing} = missing,
    args...,
) where N = begin
    object = pomp(
        object;
        rinit=rinit,
        rprocess=rprocess,
        logdmeasure=logdmeasure,
        logdprior=logdprior,
        args...,
    )
    if ismissing(estimvars)
        estimvars = keys(coef(object))
    end
    function(theta)
        traj_match_internal(theta,estimvars,object)
    end
end

traj_match_objfun(_...) = error("Incorrect call to `traj_match_objfun`.")

traj_match_internal(
    theta,
    estimvars::Union{NTuple{N,Symbol},Vector{Symbol}},
    object::AbstractPompObject,
) where N = begin
    @assert length(theta)==length(estimvars) "incorrect argument length: should be $(length(estimvars))"
    params = merge(coef(object),(;zip(estimvars,theta)...))
    x = simulate_array(object,params=params,nsim=1)
    ll = sum(logdmeasure(object,x=x,params=params))
    reg = sum(logdprior(object,params=params))
    retval = -ll-reg
    if isfinite(retval)
        retval
    else
        LogLik(Inf)
    end
end
