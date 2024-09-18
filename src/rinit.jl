export rinit

"""
    rinit(object; t0=timezero(object), params=coef(object), nsim=1)

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
rinit(
    object::PompObject;
    t0::Real = timezero(object),
    params::Union{<:NamedTuple,Vector{<:NamedTuple}} = coef(object),
    nsim::Integer = 1,
) = begin
    if isnothing(object.rinit)
        error("The *rinit* basic component is undefined.")
    end
    try
        params = val_array(params)
        [object.rinit(;params[j]...,t0=t0)
         for i ∈ 1:nsim, j ∈ eachindex(params)]
    catch e
        if isa(e,UndefKeywordError)
            error("in `rinit`: parameter " * string(e.var) * " undefined.")
        elseif hasproperty(e,:msg)
            error("in `rinit`: " * e.msg)
        else
            throw(e)            # COV_EXCL_LINE
        end
    end
end

rinit(object::AbstractPompObject;args...) = rinit(pomp(object);args...)
