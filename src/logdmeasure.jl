export logdmeasure, logdmeasure!

"""
    logdmeasure(object; times=times(object), y = obs(object), x, params)

`logdmeasure` is the workhorse for the evaluator of the log measurement density.
"""
logdmeasure(
    object::AbstractPompObject{T};
    times::Union{T,AbstractVector{T}} = times(object),
    y::Union{Y,AbstractArray{Y,M}} = obs(object),
    x::Union{X,AbstractArray{X,N}},
    params::Union{P,AbstractVector{P}},
) where {M,N,T<:Time,Y<:NamedTuple,X<:NamedTuple,P<:NamedTuple} = begin
    try
        times = val_array(times)
        params = val_array(params)
        x = val_array(x,length(params),length(times))
        y = val_array(y,length(params),length(times))
        ell = Array{Float64}(undef,size(y,1),size(x)...)
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
    logdmeasure!(object, ell; times=times(object), y = obs(object), x, params)

`logdmeasure!` is the in-place version of the `logdmeasure` workhorse.
"""
logdmeasure!(
    object::AbstractPompObject{T},
    ell::AbstractArray{Float64,4};
    times::AbstractVector{T} = times(object),
    y::AbstractArray{Y,M} = obs(object),
    x::AbstractArray{X,N},
    params::AbstractVector{P},
) where {M,N,T<:Time,Y<:NamedTuple,X<:NamedTuple,P<:NamedTuple} = begin
    try
        @assert length(params)==size(x,2)
        @assert length(params)==size(y,2)
        @assert length(times)==size(x,3)
        @assert length(times)==size(y,3)
        logdmeasure_internal!(
            pomp(object).logdmeasure,ell;
            times=times,y=y,x=x,params=params
        )
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
    ell::AbstractArray{Float64,4};
    _...,
) = begin
    for i ∈ eachindex(ell)
        @inbounds ell[i] = 0.0
    end
end

logdmeasure_internal!(
    f::Function,
    ell::AbstractArray{Float64,4};
    times::AbstractVector{T},
    y::AbstractArray{Y,3},
    x::AbstractArray{X,3},
    params::AbstractVector{P},
) where {T<:Time,Y<:NamedTuple,X<:NamedTuple,P<:NamedTuple} = begin
    for iy ∈ axes(y,1), ix ∈ axes(x,1), j ∈ eachindex(params), k ∈ eachindex(times)
        @inbounds ell[iy,ix,j,k] = f(;t=times[k],y[iy,j,k]...,x[ix,j,k]...,params[j]...)::Float64
    end
end
