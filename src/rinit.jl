"""
    rinit(object; t0=timezero(object), params=coef(object), nsim=1)

`rinit` is the workhorse for the simulator of the initial-state distribution.

## Arguments

- `object`: the AbstractPompObject
- `params`: a NamedTuple of parameters or vector of NamedTuples
- `t0`: the time at which `rinit` is to be simulated.
  This should be a single scalar.
- `nsim`: the number of simulations desired.
"""
rinit(
    object::PompObject{T,X};
    t0::T1=timezero(object),
    params::Union{P,AbstractVector{P}}=coef(object),
    nsim::Integer=1,
) where {T,X,T1<:Time,P<:NamedTuple} = begin
    params = val_array(params)
    rinit_internal(object.rinit, X, T(t0), params, object.userdata, nsim)
end

rinit(
    object::AbstractPompObject;
    kwargs...,
) = rinit(pomp(object); kwargs...)

"""
    rinit!(object, x0; t0=timezero(object), params = coef(object))

`rinit!` is the in-place version of the [`rinit`](@ref) workhorse.
"""
rinit!(
    object::PompObject{T,X},
    x0::AbstractArray{X1,2};
    t0::T1=timezero(object),
    params::Union{P,AbstractVector{P}}=coef(object),
) where {T,X,X1<:NamedTuple,T1<:Time,P<:NamedTuple} = begin
    params = val_array(params)
    rinit_internal!(x0, object.rinit, T(t0), params, object.userdata)
end

rinit!(
    object::PompObject{T,X},
    x0::AbstractArray{@NamedTuple{},2};
    t0::T1=timezero(object),
    params::Union{P,AbstractVector{P}}=coef(object),
) where {T,X,T1<:Time,P<:NamedTuple} = begin
    params = val_array(params)
    rinit_internal!(x0, object.rinit, T(t0), params, object.userdata)
end

rinit!(
    object::AbstractPompObject,
    x0::AbstractArray{X,2};
    kwargs...,
) where {X<:NamedTuple} = rinit!(pomp(object), x0; kwargs...)

rinit_internal(
    f::Nothing,
    X::Type,
    t0::Any,
    params::AbstractVector{P},
    userdata::U,
    nsim::Integer=1,
) where {P<:NamedTuple,U<:NamedTuple} =
    fill((;), length(params), nsim)

rinit_internal(
    f::Function,
    X::Type{@NamedTuple{}},
    t0::T,
    params::AbstractVector{P},
    userdata::U,
    nsim::Integer=1,
) where {T<:Time,P<:NamedTuple,U<:NamedTuple} =
    [f(; params[i]..., userdata..., t0=t0)::NamedTuple for i ∈ eachindex(params), _ ∈ 1:nsim]

rinit_internal(
    f::Function,
    X::Type,
    t0::T,
    params::AbstractVector{P},
    userdata::U,
    nsim::Integer=1,
) where {T<:Time,P<:NamedTuple,U<:NamedTuple} =
    [X(f(; params[i]..., userdata..., t0=t0)) for i ∈ eachindex(params), _ ∈ 1:nsim]

rinit_internal!(
    x0::AbstractArray{@NamedTuple{}},
    f::Nothing,
    _...,
) = nothing

rinit_internal!(
    x0::AbstractArray{X,2},
    f::Function,
    t0::T,
    params::AbstractVector{P},
    userdata::U,
) where {T<:Time,X<:NamedTuple,P<:NamedTuple,U<:NamedTuple} = begin
    for i ∈ eachindex(params), j ∈ axes(x0, 2)
        x0[i, j] = X(f(; params[i]..., userdata..., t0=t0))
    end
    nothing
end
