"""
    rprocess(object; x0, t0 = timezero(object), times=times(object), params = coef(object))

`rprocess` is the workhorse for the simulator of the process

If there is no user-supplied *rprocess* component, the dynamics are trivial.
"""
rprocess(
    object::AbstractPompObject;
    x0::Union{X,AbstractArray{X,N}} = init_state(object),
    t0::T = timezero(object),
    times::Union{T,AbstractVector{T}} = times(object),
    params::Union{P,AbstractVector{P}} = coef(object),
) where {N,T<:Time,X<:NamedTuple,P<:NamedTuple} = let
    try
        times = val_array(times)
        params = val_array(params)
        x0 = val_array(x0,length(params))
        x = similar(x0,length(times),size(x0)...)
        rprocess!(object,x,x0=x0,t0=t0,times=times,params=params)
        x
    catch e
        if hasproperty(e,:msg)
            error("in `rprocess`: " * e.msg)
        else
            throw(e)            # COV_EXCL_LINE
        end
    end
end

"""
    rprocess!(object, x; x0 = init_state(object), t0 = timezero(object), times=times(object), params = coef(object))

`rprocess!` is the in-place version of the `rprocess` workhorse.
"""
rprocess!(
    object::AbstractPompObject,
    x::AbstractArray{X,3};
    x0::AbstractArray{X,2} = init_state(object),
    t0::T = timezero(object),
    times::AbstractVector{T} = times(object),
    params::Union{P,AbstractVector{P}} = coef(object),
) where {T<:Time,X<:NamedTuple,P<:NamedTuple} = begin
    try
        params = val_array(params)
        @assert length(times)==size(x,1)
        @assert length(params)==size(x,2)
        @assert length(params)==size(x0,1)
        @assert size(x0,2)==size(x,3)
        rproc_internal!(
            x,
            pomp(object).rprocess,
            x0,times,t0,
            params,
            pomp(object).accumvars
        )
    catch e
        if isa(e,UndefKeywordError)
            error("in `rprocess!`: parameter " * string(e.var) * " undefined.")
        elseif isa(e,MethodError)
            error("in `rprocess!`: no matching method for args " * string(e.args[1]))
        elseif hasproperty(e,:msg)
            error("in `rprocess!`: " * e.msg)
        else
            throw(e)            # COV_EXCL_LINE
        end
    end
end

## the default (persistence) rprocess
rproc_internal!(
    x::AbstractArray{X,3},
    plugin::Nothing,
    x0::AbstractArray{X,2},
    times::AbstractVector{T},
    _...,
) where {T<:Time,X<:NamedTuple} = begin
    for k ∈ eachindex(times)
        @inbounds x[k,:,:] = x0
    end
end

## advance the state for each IC and parameter
## without accumulators
rproc_internal!(
    x::AbstractArray{X,3},
    plugin::PompPlugin,
    x0::AbstractArray{X,2},
    times::AbstractVector{T},
    t0::T,
    params::AbstractVector{P},
    accumvars::Nothing,
) where {T<:Time,X<:NamedTuple,P<:NamedTuple} = let
    for j ∈ eachindex(params), k ∈ axes(x0,2)
        t = t0
        @inbounds x1 = x0[j,k]
        for i ∈ eachindex(times)
            @inbounds t,x1 = rprocess_step(plugin,t,times[i],x1,params[j])
            @inbounds x[i,j,k] = x1
        end
    end
end

## advance the state for each IC and parameter
## without accumulators
rproc_internal!(
    x::AbstractArray{X,3},
    plugin::PompPlugin,
    x0::AbstractArray{X,2},
    times::AbstractVector{T},
    t0::T,
    params::AbstractVector{P},
    accumvars::A
) where {T<:Time,X<:NamedTuple,P<:NamedTuple,A<:NamedTuple} = let
    for j ∈ eachindex(params), k ∈ axes(x0,2)
        t = t0
        @inbounds x1 = x0[j,k]
        for i ∈ eachindex(times)
            x1 = merge(x1,accumvars)::X
            @inbounds t,x1 = rprocess_step(plugin,t,times[i],x1,params[j])
            @inbounds x[i,j,k] = x1
        end
    end
end
