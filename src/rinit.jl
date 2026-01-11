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
) where {T,P<:NamedTuple} = begin
    params = val_array(params)
    rinit_internal(pomp(object).rinit,t0,params,nsim)
end

"""
    rinit!(object, x0; t0=timezero(object), params = coef(object))

`rinit!` is the in-place version of the `rinit` workhorse.
"""
rinit!(
    object::AbstractPompObject,
    x0::AbstractArray{X,2};
    t0::T = timezero(object),
    params::Union{P,AbstractVector{P}} = coef(object),
) where {T,X,P<:NamedTuple} = begin
    params = val_array(params)
    rinit_internal!(x0,pomp(object).rinit,t0,params)
end

rinit_internal(
    f::Nothing,
    t0::Any,
    params::AbstractVector{P},
    nsim::Integer = 1
) where {P<:NamedTuple} =
    fill((;),length(params),nsim)

rinit_internal(
    f::Function,
    t0::T,
    params::AbstractVector{P},
    nsim::Integer = 1,
) where {T<:Time,P<:NamedTuple} =
    [f(;params[i]...,t0=t0) for i ∈ eachindex(params), _ ∈ 1:nsim]

rinit_internal!(                # COV_EXCL_LINE
    x0::AbstractArray{X},
    f::Nothing,
    _...,
) where {X} = begin
    fill!(x0,(;))
    nothing
end

rinit_internal!(
    x0::AbstractArray{X,2},
    f::Function,
    t0::T,
    params::AbstractVector{P},
) where {T<:Time,X,P<:NamedTuple} = begin
    for i ∈ eachindex(params), j ∈ axes(x0,2)
        x0[i,j] = f(;params[i]...,t0=t0)::X
    end                         # COV_EXCL_LINE
end
