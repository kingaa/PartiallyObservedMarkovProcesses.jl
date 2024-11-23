export rprocess, rprocess!

"""
    rprocess(object; x0, t0 = timezero(object), times=times(object), params)

`rprocess` is the workhorse for the simulator of the process

If there is no user-supplied *rprocess* component, the dynamics are trivial.
"""
rprocess(
    object::AbstractPompObject{T};
    x0::Union{X,AbstractArray{X,N}},
    t0::T = timezero(object),
    times::Union{T,AbstractVector{T}} = times(object),
    params::Union{P,AbstractVector{P}},
) where {N,T<:Time,X<:NamedTuple,P<:NamedTuple} = let
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

"""
    rprocess!(object, x; x0, t0 = timezero(object), times=times(object), params)

`rprocess!` is the in-place version of the `rprocess` workhorse.
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
    f::Nothing,
    x0::AbstractArray{X,2},
    times::AbstractVector{T},
    _...,
) where {T<:Time,X<:NamedTuple} = begin
    for k ∈ eachindex(times)
        @views x[:,:,k] = x0
    end
end

## advance the state for each IC and parameter
## this the case with no accumulator variables
rproc_internal!(
    x::AbstractArray{X,3},
    f::Function,
    x0::AbstractArray{X,2},
    times::AbstractVector{T},
    t0::T,
    params::AbstractVector{P},
    accumvars::Nothing,
) where {T<:Time,X<:NamedTuple,P<:NamedTuple} = let
    for i ∈ axes(x0,1), j ∈ eachindex(params)
        t = t0
        x1 = x0[i,j]
        for k ∈ eachindex(times)
            while t < times[k]
                t,x1... = f(;t=t,x1...,params[j]...)
            end
            x[i,j,k] = x1
        end
    end
end

## advance the state for each IC and parameter
## this the case with accumulator variables
rproc_internal!(
    x::AbstractArray{X,3},
    f::Function,
    x0::AbstractArray{X,2},
    times::AbstractVector{T},
    t0::T,
    params::AbstractVector{P},
    accumvars::NamedTuple{N},
) where {M,N,T<:Time,X<:NamedTuple{M},P<:NamedTuple} = let
    local Q = Tuple(setdiff(M,N)) # non-accumulator variables
    for i ∈ axes(x0,1), j ∈ eachindex(params)
        t = t0
        x1 = x0[i,j]
        for k ∈ eachindex(times)
            x1 = merge(x1,accumvars)
            while t < times[k]
                t,x1... = f(;t=t,x1...,params[j]...)
            end
            x[i,j,k] = x1
        end
    end
end
