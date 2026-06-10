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

Note that most errors that arise in the integration of the vectorfield will be
trapped, with a warning. The objective function takes values Inf in this case.
"""
traj_match_objfun(
    object::ValidPompData,
    estimvars::Union{NTuple{N,Symbol},Missing} = missing;
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
    estimvars::NTuple{N,Symbol},
    object::AbstractPompObject,
) where N = begin
    @assert length(theta)==N "incorrect argument length: should be $N"
    params = merge(coef(object),(;zip(estimvars,theta)...))
    try
        x = simulate_array(object,params=params,nsim=1)
        ll = sum(logdmeasure(object,x=x,params=params))
        reg = sum(logdprior(object,params=params))
        retval = -ll-reg
        if isfinite(retval)
            retval
        else
            LogLik(Inf)
        end
    catch e
        if e isa UndefKeywordError
            throw(e)
        elseif hasproperty(e,:msg)
            @warn("in `traj_match_objfun`: $(e.msg)")
        else
            @warn("in `traj_match_objfun`: error: $e")
        end
        LogLik(Inf)
    end
end
