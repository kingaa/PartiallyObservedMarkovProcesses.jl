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
    rinit_internal(object.rinit, T(t0), params, object.userdata, nsim)
end

rinit(
    object::AbstractPompObject;
    args...,
) = rinit(pomp(object);args...)

"""
    rinit!(object, x0; t0=timezero(object), params = coef(object))

`rinit!` is the in-place version of the [`rinit`](@ref) workhorse.
"""
rinit!(
    object::PompObject{T,X},
    x0::AbstractArray{X1,2};
    t0::T1=timezero(object),
    params::Union{P,AbstractVector{P}}=coef(object),
) where {T,X,T1<:Time,X1<:NamedTuple,P<:NamedTuple} = begin
    params = val_array(params)
    rinit_internal!(x0, object.rinit, T(t0), params, object.userdata)
end

rinit!(
    object::AbstractPompObject,
    x0::AbstractArray{X,2};
    args...,
) where {X<:NamedTuple} = rinit!(pomp(object),x0;args...)

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

# COV_EXCL_START  (to bypass bug in LocalCoverage.jl)
rinit_internal!(
    x0::AbstractArray{X},
    f::Nothing,
    _...,
) where {X<:NamedTuple} = begin
    # COV_EXCL_STOP
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
