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
    t0::T=timezero(object),
    params::Union{P,AbstractVector{P}}=coef(object),
    nsim::Integer=1,
) where {T<:Time,P<:NamedTuple} = begin
    params = val_array(params)
    rinit_internal(pomp(object).rinit, t0, params, pomp(object).userdata, nsim)
end

"""
    rinit!(object, x0; t0=timezero(object), params = coef(object))

`rinit!` is the in-place version of the `rinit` workhorse.
"""
rinit!(
    object::AbstractPompObject,
    x0::AbstractArray{X,2};
    t0::T=timezero(object),
    params::Union{P,AbstractVector{P}}=coef(object),
) where {T<:Time,X<:NamedTuple,P<:NamedTuple} = begin
    params = val_array(params)
    rinit_internal!(x0, pomp(object).rinit, t0, params, pomp(object).userdata)
end

rinit_internal(
    f::Nothing,
    t0::Any,
    params::AbstractVector{P},
    userdata::U,
    nsim::Integer=1,
) where {P<:NamedTuple,U<:NamedTuple} =
    fill((;), length(params), nsim)

rinit_internal(
    f::Function,
    t0::T,
    params::AbstractVector{P},
    userdata::U,
    nsim::Integer=1,
) where {T<:Time,P<:NamedTuple,U<:NamedTuple} =
    [f(; params[i]..., userdata..., t0=t0)::NamedTuple for i ∈ eachindex(params), _ ∈ 1:nsim]

rinit_internal!(                # COV_EXCL_LINE
    x0::AbstractArray{X},
    f::Nothing,
    _...,
) where {X<:NamedTuple} = begin
    fill!(x0, (;))
    nothing
end

rinit_internal!(
    x0::AbstractArray{X,2},
    f::Function,
    t0::T,
    params::AbstractVector{P},
    userdata::U,
) where {T<:Time,X<:NamedTuple,P<:NamedTuple,U<:NamedTuple} = begin
    for i ∈ eachindex(params), j ∈ axes(x0, 2)
        x0[i, j] = X(f(; params[i]..., userdata..., t0=t0))
    end                         # COV_EXCL_LINE
end
