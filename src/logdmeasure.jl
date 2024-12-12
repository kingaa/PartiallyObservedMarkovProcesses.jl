export logdmeasure, logdmeasure!

"""
    logdmeasure(object; times=times(object), y=obs(object), x=states(object), params=coef(object))

`logdmeasure` is the workhorse for the evaluator of the log measurement density.
"""
logdmeasure(
    object::AbstractPompObject;
    times::Union{T,AbstractVector{T}} = times(object),
    y::Union{Y,AbstractArray{Y}} = obs(object),
    x::Union{X,AbstractArray{X}} = states(object),
    params::Union{P,AbstractVector{P}} = coef(object),
) where {T<:Time,Y<:NamedTuple,X<:NamedTuple,P<:NamedTuple} = begin
    try
        times = val_array(times)
        params = val_array(params)
        x = val_array(x,length(times),length(params))
        y = val_array(y,length(times),length(params))
        ell = Array{LogLik}(undef,size(x)...,size(y,3))
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
    logdmeasure!(object, ell; times=times(object), y=obs(object), x=states(object), params=coef(object))

`logdmeasure!` is the in-place version of the `logdmeasure` workhorse.
"""
logdmeasure!(
    object::AbstractPompObject,
    ell::AbstractArray{W,4};
    times::AbstractVector{T} = times(object),
    y::AbstractArray{Y} = obs(object),
    x::AbstractArray{X} = states(object),
    params::Union{P,AbstractVector{P}} = coef(object),
) where {W<:Real,T<:Time,Y<:NamedTuple,X<:NamedTuple,P<:NamedTuple} = begin
    try
        params = val_array(params)
        @assert length(times)==size(x,1)
        @assert length(times)==size(y,1)
        @assert length(params)==size(x,2)
        @assert length(params)==size(y,2)
        logdmeasure_internal!(ell,pomp(object).logdmeasure,times,y,x,params)
    catch e
        if isa(e,UndefKeywordError)
            error("in `logdmeasure!`: parameter " * string(e.var) * " undefined.")
        elseif isa(e,MethodError)
            error("in `logdmeasure!`: no matching method for args " * string(e.args[1]))
        elseif hasproperty(e,:msg)
            error("in `logdmeasure!`: " * e.msg)
        else
            throw(e)            # COV_EXCL_LINE
        end
    end
end

logdmeasure_internal!(
    ell::AbstractArray{W,4},
    f::Nothing,
    _...,
) where {W<:Real} = begin
    for i ∈ eachindex(ell)
        @inbounds ell[i] = W(0)
    end
end

logdmeasure_internal!(
    ell::AbstractArray{W,4},
    f::Function,
    times::AbstractVector{T},
    y::AbstractArray{Y,3},
    x::AbstractArray{X,3},
    params::AbstractVector{P},
) where {W<:Real,T<:Time,Y<:NamedTuple,X<:NamedTuple,P<:NamedTuple} = begin
    for i ∈ eachindex(times), j ∈ eachindex(params), kx ∈ axes(x,3), ky ∈ axes(y,3)
        @inbounds ell[i,j,kx,ky] = f(;t=times[i],y[i,j,ky]...,x[i,j,kx]...,params[j]...)::W
    end
end
