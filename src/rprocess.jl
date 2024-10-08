export rprocess, rprocess!

rproc_internal(
    f::Nothing,
    x::AbstractArray{X,3};
    x0::AbstractArray{X,2},
    times::AbstractVector{T} = times(object),
    _...,
) where {T<:Time,X<:NamedTuple} = begin
    for k ∈ eachindex(times)
        x[:,:,k] = x0
    end
end

rproc_internal(
    f::Function,
    x::AbstractArray{X,3};
    x0::AbstractArray{X,2},
    t0::T = timezero(object),
    times::AbstractVector{T} = times(object),
    params::AbstractVector{P},
    accumvars::NamedTuple{N},
) where {N,T<:Time,X<:NamedTuple,P<:NamedTuple} = begin
    regvar = Tuple(setdiff(keys(x0[begin]),N))
    for i ∈ axes(x0,1), j ∈ eachindex(params)
        t::T = t0
        x1::X = x0[i,j]
        for k ∈ eachindex(times)
            t,x1 = advance_rproc(x1,f,t,times[k],params[j],Val(accumvars),Val(regvar))
            x[i,j,k] = x1
        end
    end
end

advance_rproc(
    x::NamedTuple{M},
    f::Function,
    t::T,
    tf::T,
    p::P,
    ::Val{Z},
    ::Val{Q},
) where {T<:Time,P<:NamedTuple,M,Z,Q} = begin
    x1 = (;x[Q]...,Z...)
    while t < tf
        (t,x1...) = f(;t=t,x1...,p...)
    end
    t,x1
end

"""
    rprocess(object; x0, t0 = timezero(object), times=times(object), params=coef(object))

`rprocess` is the workhorse for the simulator of the process

The user can supply an *rprocess* component as a function that takes states, parameters, and current time (`t`) and returns the updated time and state.

If there is no user-supplied *rprocess* component, the dynamics are trivial.
"""
rprocess!(
    object::AbstractPompObject{T},
    x::AbstractArray{X,3};
    x0::AbstractArray{X,2},
    t0::T = timezero(object),
    times::AbstractVector{T} = times(object),
    params::AbstractVector{P},
) where {T<:Time,X<:NamedTuple,P<:NamedTuple} = begin
    try
        @assert size(x0,1)==size(x,1)
        @assert length(params)==size(x0,2)
        @assert length(params)==size(x,2)
        @assert length(times)==size(x,3)
        rproc_internal(
            pomp(object).rprocess,x;
            x0=x0,t0=t0,times=times,
            params=params,
            accumvars=pomp(object).accumvars
        )
    catch e
        if isa(e,UndefKeywordError)
            error("in `rprocess!`: parameter " * string(e.var) * " undefined.")
        elseif hasproperty(e,:msg)
            error("in `rprocess!`: " * e.msg)
        else
            throw(e)            # COV_EXCL_LINE
        end
    end
end

rprocess(
    object::AbstractPompObject{T};
    x0::Union{X,AbstractArray{X,N}},
    t0::T = timezero(object),
    times::Union{T,AbstractVector{T}} = times(object),
    params::Union{P,AbstractVector{P}},
) where {N,T<:Time,X<:NamedTuple,P<:NamedTuple} = begin
    try
        times = val_array(times)
        params = val_array(params)
        x0 = val_array(x0,length(params))
        x = Array{X}(undef,size(x0)...,length(times))
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
