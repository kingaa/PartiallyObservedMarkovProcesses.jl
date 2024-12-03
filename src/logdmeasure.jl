export logdmeasure, logdmeasure!

"""
    logdmeasure(object; times=times(object), y=obs(object), x=states(object), params=coef(object))

`logdmeasure` is the workhorse for the evaluator of the log measurement density.
"""
logdmeasure(
    object::AbstractPompObject;
    times::Union{T,AbstractVector{T}} = times(object),
    y::Union{Y,AbstractArray{Y,M}} = obs(object),
    x::Union{X,AbstractArray{X,N}} = states(object),
    params::Union{P,AbstractVector{P}} = coef(object),
) where {M,N,T<:Time,Y<:NamedTuple,X<:NamedTuple,P<:NamedTuple} = begin
    try
        times = val_array(times)
        params = val_array(params)
        x = val_array(x,length(params),length(times))
        y = val_array(y,length(params),length(times))
        ell = Array{LogLik}(undef,size(y,1),size(x)...)
        logdmeasure!(object,ell;times=times,y=y,x=x,params=params)
        ell
    catch e
        if hasproperty(e,:msg)
            error("in `logdmeasure`: " * e.msg)
        else
            throw(e)            # COV_EXCL_LINE
        end
    end
end

"""
    logdmeasure!(object; times=times(object), y=obs(object), x=states(object), params=coef(object))

`logdmeasure!` is the in-place version of the `logdmeasure` workhorse.
"""
logdmeasure!(
    object::AbstractPompObject,
    ell::AbstractArray{W,4};
    times::AbstractVector{T} = times(object),
    y::AbstractArray{Y,M} = obs(object),
    x::AbstractArray{X,N} = states(object),
    params::Union{P,AbstractVector{P}} = coef(object),
) where {M,N,W<:Real,T<:Time,Y<:NamedTuple,X<:NamedTuple,P<:NamedTuple} = begin
    try
        params = val_array(params)
        @assert length(params)==size(x,2)
        @assert length(params)==size(y,2)
        @assert length(times)==size(x,3)
        @assert length(times)==size(y,3)
        logdmeasure_internal!(pomp(object).logdmeasure,ell,times,y,x,params)
    catch e
        if isa(e,UndefKeywordError)
            error("in `logdmeasure!`: parameter " * string(e.var) * " undefined.")
        elseif hasproperty(e,:msg)
            error("in `logdmeasure!`: " * e.msg)
        else
            throw(e)            # COV_EXCL_LINE
        end
    end
end

logdmeasure_internal!(
    f::Nothing,
    ell::AbstractArray{W,4},
    _...,
) where {W<:Real} = begin
    for i ∈ eachindex(ell)
        @inbounds ell[i] = W(0.0)
    end
end

logdmeasure_internal!(
    f::Function,
    ell::AbstractArray{W,4},
    times::AbstractVector{T},
    y::AbstractArray{Y,3},
    x::AbstractArray{X,3},
    params::AbstractVector{P},
) where {W<:Real,T<:Time,Y<:NamedTuple,X<:NamedTuple,P<:NamedTuple} = begin
    for iy ∈ axes(y,1), ix ∈ axes(x,1), j ∈ eachindex(params), k ∈ eachindex(times)
        @inbounds ell[iy,ix,j,k] = f(;t=times[k],y[iy,j,k]...,x[ix,j,k]...,params[j]...)::W
    end
end
