"""
    rprior(object; params=coef(object), nsim = 1)

`rprior` is the workhorse for the simulator of the prior distribution.
If `nsim > 1`, then a matrix is returned.
"""
rprior(
    object::AbstractPompObject;
    params::Union{P,AbstractVector{P}} = coef(object),
    nsim::Integer = 1,
) where {P<:NamedTuple} = begin
    params = val_array(params)
    rprior_internal(pomp(object).rprior,params,nsim,pomp(object).userdata)
end

rprior_internal(
    f::Nothing,
    params::AbstractVector{P},
    nsim::Integer,
    _...,
) where {P<:NamedTuple} = begin
    @assert nsim > 0
    repeat(params,1,nsim)
end

rprior_internal(
    f::Function,
    params::AbstractVector{P},
    nsim::Integer,
    userdata::U,
) where {P<:NamedTuple,U<:NamedTuple} = begin
    @inbounds(
        [f(;params[j]...,userdata...)::NamedTuple
         for j ∈ eachindex(params), _ ∈ 1:nsim]
    )
end
