export rinit

rinit_internal(
    f::Nothing;
    params::AbstractVector{P},
    nsim::Integer = 1,
    _...,
) where {P<:NamedTuple} = begin
    reshape(fill((),nsim*length(params)),nsim,length(params))
end

rinit_internal(
    f::Function;
    t0::T,
    params::AbstractVector{P},
    nsim::Integer = 1,
) where {T,P<:NamedTuple} = begin
    [f(;params[j]...,t0=t0)
     for i ∈ 1:nsim, j ∈ eachindex(params)]
end


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
"""
rinit(
    object::AbstractPompObject{T};
    t0::T = timezero(object),
    params::Union{P,AbstractVector{P}},
    nsim::Integer = 1,
) where {T,P<:NamedTuple} = begin
    try
        params = val_array(params)
        rinit_internal(pomp(object).rinit,t0=t0,params=params,nsim=nsim)
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
