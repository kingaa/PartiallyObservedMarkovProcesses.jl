export rinit

"""
    rinit(object, params=object.params, t0=object.t0, nsim=1)

`rinit` is the workhorse for the simulator of the initial-state distribution.

The user can supply an *rinit* component as a function that takes parameters and, optionally, `t0`, the initial time.

## Arguments

- `object`: the PompObject
- `params`: a NamedTuple of parameters
- `t0`: the time at which `rinit` is to be simulated.
  This should be a single scalar.
- `nsim`: the number of simulations desired.

Calling `rinit()` in the absence of a user-supplied *rinit* component results in an error.
"""
rinit = function (
    object::PompObject;
    params::NamedTuple = object.params,
    t0::Real = object.t0,
    nsim::Integer = 1
    )
    if isnothing(object.rinit)
        error("The *rinit* basic component is undefined.")
    end
    try
        [object.rinit(;params...,t0=t0) for _ ∈ 1:nsim]
    catch e
        if isa(e,UndefKeywordError)
            error("in `rinit`: parameter " * e.var * " undefined.")
        elseif hasproperty(e,:msg)
            error("in `rinit`: " * e.msg)
        else
            throw(e)
        end
    end
end
