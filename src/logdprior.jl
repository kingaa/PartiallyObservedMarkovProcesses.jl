"""
    logdprior(object; params=coef(object))

`logdprior` is the workhorse for the evaluator of the log prior density.
If no prior is specified, the `logdprior` returns 0 for all inputs.
"""
logdprior(
    object::AbstractPompObject;
    params::Union{P,AbstractVector{P}}=coef(object),
) where {P<:NamedTuple} = begin
    params = val_array(params)
    ell = Array{LogLik}(undef, size(params))
    logdprior!(object, ell; params=params)
    ell
end

"""
    logdprior!(object, ell; params=coef(object))

`logdprior!` is the in-place version of the `logdprior` workhorse.
"""
logdprior!(
    object::AbstractPompObject,
    ell::AbstractArray{W,1};
    params::Union{P,AbstractVector{P}}=coef(object),
) where {W<:AbstractFloat,P<:NamedTuple} = begin
    params = val_array(params)
    logdprior_internal!(ell, pomp(object).logdprior, params, pomp(object).userdata)
end

logdprior_internal!(
    ell::AbstractArray{W,1},
    f::Nothing,
    _...,
) where {W<:AbstractFloat} = begin
    for i ∈ eachindex(ell)
        @inbounds ell[i] = W(0)
    end
end

logdprior_internal!(
    ell::AbstractArray{W,1},
    f::Function,
    params::AbstractVector{P},
    userdata::U,
) where {W<:AbstractFloat,P<:NamedTuple,U<:NamedTuple} = begin
    for j ∈ eachindex(params)
        @inbounds ell[j] = W(f(; params[j]..., userdata...))
    end
end
