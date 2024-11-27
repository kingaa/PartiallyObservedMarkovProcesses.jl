export rinit, rinit!

"""
    rinit(object; t0=timezero(object), params=coef(object), nsim=1)

`rinit` is the workhorse for the simulator of the initial-state distribution.

## Arguments

- `object`: the PompObject
- `params`: a NamedTuple of parameters or vector of NamedTuples
- `t0`: the time at which `rinit` is to be simulated.
  This should be a single scalar.
- `nsim`: the number of simulations desired.
"""
rinit(
    object::AbstractPompObject;
    t0::T = timezero(object),
    params::Union{P,AbstractVector{P}} = coef(object),
    nsim::Integer = 1,
) where {T,P<:NamedTuple} = let
    try
        params = val_array(params)
        rinit_internal(pomp(object).rinit,t0,params,nsim)
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

"""
    rinit!(object, x0; t0=timezero(object), params)

`rinit!` is the in-place version of the `rinit` workhorse.
"""
rinit!(
    object::AbstractPompObject,
    x0::AbstractArray{X,2};
    t0::T = timezero(object),
    params::Union{P,AbstractVector{P}} = coef(object),
) where {T,X,P<:NamedTuple} = let
    try
        params = val_array(params)
        rinit_internal!(pomp(object).rinit,x0,t0,params)
    catch e
        if isa(e,UndefKeywordError)
            error("in `rinit!`: parameter " * string(e.var) * " undefined.")
        elseif hasproperty(e,:msg)
            error("in `rinit!`: " * e.msg)
        else
            throw(e)            # COV_EXCL_LINE
        end
    end
end

rinit_internal(
    f::Nothing,
    t0::Any,
    params::AbstractVector{P},
    nsim::Integer = 1
) where {P<:NamedTuple} =
    fill((;),nsim,length(params))

rinit_internal(
    f::Function,
    t0::T,
    params::AbstractVector{P},
    nsim::Integer = 1,
) where {T<:Time,P<:NamedTuple} =
    [f(;params[j]...,t0=t0) for i ∈ 1:nsim, j ∈ eachindex(params)]

rinit_internal!(
    f::Nothing,
    x0::AbstractArray{X},
    t0::Any,
    params::AbstractVector{P},
) where {X,P<:NamedTuple} = begin
    fill!(x0,(;))
    nothing
end

rinit_internal!(
    f::Function,
    x0::AbstractArray{X,2},
    t0::T,
    params::AbstractVector{P},
) where {T<:Time,X,P<:NamedTuple} = begin
    for i ∈ axes(x0,1), j ∈ eachindex(params)
        x0[i,j] = f(;params[j]...,t0=t0)::X
    end                         # COV_EXCL_LINE
    nothing                     # COV_EXCL_LINE
end
