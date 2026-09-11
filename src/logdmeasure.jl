"""
    logdmeasure(object; times=times(object), y=obs(object),
                x=states(object), params=coef(object))

`logdmeasure` is the workhorse for the evaluator of the log
measurement density.
"""
logdmeasure(
    object::AbstractPompObject;
    times::Union{T,AbstractVector{T}}=times(object),
    y::Union{Y,AbstractVector{Y}}=obs(object),
    x::Union{X,AbstractArray{X}}=states(object),
    params::Union{P,AbstractVector{P}}=coef(object),
) where {T<:Time,Y<:NamedTuple,X<:NamedTuple,P<:NamedTuple} = begin
    times = val_array(times)
    y = val_array(y)
    params = val_array(params)
    x = val_array(x, length(times), length(params))
    @assert length(y)==length(times)
    ell = Array{LogLik}(undef, size(x))
    logdmeasure_internal!(
        ell,
        pomp(object).logdmeasure,
        times,
        y, x, params,
        pomp(object).userdata
    )
    ell
end

"""
    logdmeasure!(object, ell; times=times(object), y=obs(object),
                 x=states(object), params=coef(object))

`logdmeasure!` is the in-place version of the `logdmeasure` workhorse.
If no `logdmeasure` component has been specified, this returns 0 for
all inputs.
"""
logdmeasure!(
    object::AbstractPompObject,
    ell::AbstractArray{W,3};
    times::Union{T,AbstractVector{T}}=times(object),
    y::Union{Y,AbstractVector{Y}}=obs(object),
    x::Union{X,AbstractArray{X}}=states(object),
    params::Union{P,AbstractVector{P}}=coef(object),
) where {
    W<:AbstractFloat,T<:Time,Y<:NamedTuple,
    X<:NamedTuple,P<:NamedTuple
} = begin
    times = val_array(times)
    y = val_array(y)
    params = val_array(params)
    x = val_array(x, length(times), length(params))
    @assert length(y)==length(times)
    @assert size(ell)==size(x)
    logdmeasure_internal!(
        ell,
        pomp(object).logdmeasure,
        times,
        y, x, params,
        pomp(object).userdata
    )
    nothing
end

# COV_EXCL_START  (to bypass bug in LocalCoverage.jl)
logdmeasure_internal!(
    ell::AbstractArray{W,3},
    f::Nothing,
    _...,
) where {W<:AbstractFloat} = begin
    # COV_EXCL_STOP
    for i ∈ eachindex(ell)
        ell[i] = W(0)
    end
    nothing
end

logdmeasure_internal!(
    ell::AbstractArray{W,3},
    f::Function,
    times::AbstractArray{T,1},
    y::AbstractArray{Y,1},
    x::AbstractArray{X,3},
    params::AbstractArray{P,1},
    userdata::U,
) where {
    W<:AbstractFloat,T<:Time,Y<:NamedTuple,X<:NamedTuple,
    P<:NamedTuple,U<:NamedTuple
} = begin
    for i ∈ axes(ell,1), j ∈ axes(ell,2), k ∈ axes(ell,3)
        ell[i,j,k] = W(f(; t=times[i], y[i]..., x[i,j,k]..., params[j]..., userdata...))
    end
    nothing
end
