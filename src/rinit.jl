export rinit, rinit!

"""
    rinit(object; t0=timezero(object), params, nsim=1)

`rinit` is the workhorse for the simulator of the initial-state distribution.

## Arguments

- `object`: the PompObject
- `params`: a NamedTuple of parameters
- `t0`: the time at which `rinit` is to be simulated.
  This should be a single scalar.
- `nsim`: the number of simulations desired.
"""
rinit(
    object::AbstractPompObject{T};
    t0::T = timezero(object),
    params::Union{P,AbstractVector{P}},
    nsim::Integer = 1,
) where {T,P<:NamedTuple} = begin
    try
        params = val_array(params)
        rinit_internal(pomp(object).rinit,t0=t0,params=params,nsim=nsim)
    catch e
        if isa(e,UndefKeywordError)
            error("in `rinit`: parameter " * string(e.var) * " undefined.")
        elseif hasproperty(e,:msg)
            error("in `rinit`: " * e.msg)
        else
            throw(e)            # COV_EXCL_LINE
        end
    end
end

"""
    rinit!(object, x0; t0=timezero(object), params)

`rinit!` is the in-place version of the `rinit` workhorse.
"""
rinit!(
    object::AbstractPompObject{T},
    x0::AbstractArray{X,2};
    t0::T = timezero(object),
    params::Union{P,AbstractVector{P}},
) where {T,X,P<:NamedTuple} = begin
    try
        params = val_array(params)
        rinit_internal!(pomp(object).rinit,x0,t0=t0,params=params)
    catch e
        if isa(e,UndefKeywordError)
            error("in `rinit!`: parameter " * string(e.var) * " undefined.")
        elseif hasproperty(e,:msg)
            error("in `rinit!`: " * e.msg)
        else
            throw(e)            # COV_EXCL_LINE
        end
    end
end

rinit_internal(
    f::Nothing;
    params::AbstractVector{P},
    nsim::Integer = 1,
    _...,
) where {P<:NamedTuple} = begin
    fill((;),nsim,length(params))
end

rinit_internal(
    f::Function;
    t0::T,
    params::AbstractVector{P},
    nsim::Integer = 1,
) where {T,P<:NamedTuple} = begin
    [f(;params[j]...,t0=t0)
     for i ∈ 1:nsim, j ∈ eachindex(params)]
end

rinit_internal!(
    f::Nothing,
    x0::AbstractArray{X};
    params::AbstractVector{P},
    _...,
) where {X,P<:NamedTuple} = begin
    fill!(x0,(;))
end

rinit_internal!(
    f::Function,
    x0::AbstractArray{X,2};
    t0::T,
    params::AbstractVector{P},
) where {T,X,P<:NamedTuple} = begin
    for i ∈ axes(x0,1), j ∈ eachindex(params)
        x0[i,j] = f(;params[j]...,t0=t0)::X
    end
end
