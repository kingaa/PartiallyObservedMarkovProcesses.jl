"""
    traj_match_objfun(object, estimvars; rinit, rprocess,
        logdmeasure, logdprior, bigvalue = Inf, whitelist,
        kwargs...)

Returns a function that compares a model trajectory with the data,
returning minus the sum of the log likelihood and log prior density.
Trajectory matching consists of the minimization of such an objective
function.  At least the `rinit`, `rprocess`, and `logdmeasure` basic
components are needed. `kwargs...` can be used to modify or unset
these fields.  `rprocess` must be a `VectorfieldPlugin` (i.e.,
constructed via a call to [`vectorfield`](@ref)).

The returned function takes as input a vector or tuple of length equal
to that of `estimvars`. This vector is associated,
element-for-element, with the symbols in `estimvars`. By default,
`estimvars = keys(coef(object))`.

Note that failure of the integrator (signaled via an unsuccessful
`retcode`) will result in an objective-function value of `bigvalue`
and a warning message. Depending on the nature of the objective
function, and on the integration algorithm, additional exceptions may
arise. In some cases, one may wish to trap these and return an
objective-function value of `bigvalue`. The `whitelist` argument
specifies an exception type (or `Union` of exception types) that
should be handled in this way.
"""
traj_match_objfun(
    object::ValidPompData,
    estimvars::Union{NTuple{N,Symbol},Vector{Symbol},Missing} = missing;
    rinit::Union{Function,Nothing,Missing} = missing,
    rprocess::Union{VectorfieldPlugin,Nothing,Missing} = missing,
    logdmeasure::Union{Function,Nothing,Missing} = missing,
    logdprior::Union{Function,Nothing,Missing} = missing,
    whitelist::Type = Nothing,
    bigvalue::AbstractFloat = Inf,
    kwargs...,
) where N = begin
    object = pomp(
        object;
        rinit=rinit,
        rprocess=rprocess,
        logdmeasure=logdmeasure,
        logdprior=logdprior,
        kwargs...,
    )
    if ismissing(estimvars)
        estimvars = keys(coef(object))
    end
    function(theta)
        traj_match_internal(theta,estimvars,object,LogLik(bigvalue),whitelist,)
    end
end

traj_match_objfun(_...) = error("Incorrect call to `traj_match_objfun`.")

traj_match_internal(theta, estimvars, object, bigvalue, whitelist,) = begin
    @assert length(theta)==length(estimvars) "incorrect argument length: should be $(length(estimvars))"
    params = merge(coef(object),(;zip(estimvars,theta)...))
    try
        x = simulate_array(object,params=params,nsim=1)
        ll = sum(logdmeasure(object,x=x,params=params))
        reg = sum(logdprior(object,params=params))
        retval = -ll-reg
        if isfinite(retval)
            retval
        else
            bigvalue
        end
    catch e
        if e isa FailedIntegrationException
            @warn("in trajectory matching: $e")
        elseif e isa whitelist
            @error("in trajectory matching: whitelisted error: $e")
        else
            throw(e)
        end
        bigvalue
    end
end
