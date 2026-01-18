"""
    logdprior(object; params=coef(object))

`logdprior` is the workhorse for the evaluator of the log prior density.
"""
logdprior(
    object::AbstractPompObject;
    params::Union{P,AbstractVector{P}} = coef(object),
) where {P<:NamedTuple} = begin
    params = val_array(params)
    ell = Array{LogLik}(undef,size(params))
    logdprior!(object,ell;params=params)
    ell
end

"""
    logdprior!(object, ell; params=coef(object))

`logdprior!` is the in-place version of the `logdprior` workhorse.
"""
logdprior!(
    object::AbstractPompObject,
    ell::AbstractArray{W,1};
    params::Union{P,AbstractVector{P}} = coef(object),
) where {W<:Real,P<:NamedTuple} = begin
    params = val_array(params)
    logdprior_internal!(ell,pomp(object).logdprior,params,pomp(object).userdata)
end

logdprior_internal!(
    ell::AbstractArray{W,1},
    f::Nothing,
    _...,
) where {W<:Real} = begin
    for i ∈ eachindex(ell)
        @inbounds ell[i] = W(0)
    end
end

logdprior_internal!(
    ell::AbstractArray{W,1},
    f::Function,
    params::AbstractVector{P},
    userdata::U,
) where {W<:Real,P<:NamedTuple,U<:NamedTuple} = begin
    for j ∈ eachindex(params)
        @inbounds ell[j] = f(;params[j]...,userdata...)::W
    end
end
