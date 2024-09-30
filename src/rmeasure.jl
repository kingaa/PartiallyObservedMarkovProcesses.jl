export rmeasure

rmeas_internal(
    f::Nothing,
    x::AbstractArray{X,3};
    _...,
) where {X<:NamedTuple} = begin
    reshape(fill((),length(x)),size(x)...)
end

rmeas_internal(
    f::Function,
    x::AbstractArray{X,3};
    times::AbstractVector{T},
    params::AbstractVector{P},
) where {T,X<:NamedTuple,P<:NamedTuple} = begin
    [f(;t=times[k],x[i,j,k]...,params[j]...)
     for i ∈ axes(x,1), j ∈ eachindex(params), k ∈ eachindex(times)]
end


"""
    rmeasure(object; x, times=times(object), params=coef(object))

`rmeasure` is the workhorse for the simulator of the measurement distribution.

The user can supply an *rmeasure* component as a function that takes states, parameters, and, optionally, `t`, the current time.
"""
rmeasure(
    object::AbstractPompObject{T};
    x::Union{X,AbstractArray{X,N}},
    times::Union{T,AbstractVector{T}} = times(object),
    params::Union{P,AbstractVector{P}},
) where {N,T,X<:NamedTuple,P<:NamedTuple} = begin
    try
        times = val_array(times)
        params = val_array(params)
        x = val_array(x,length(params),length(times))
        rmeas_internal(pomp(object).rmeasure,x,times=times,params=params)
    catch e
        if isa(e,UndefKeywordError)
            error("in `rmeasure`: parameter " * string(e.var) * " undefined.")
        elseif hasproperty(e,:msg)
            error("in `rmeasure`: " * e.msg)
        else
            throw(e)            # COV_EXCL_LINE
        end
    end
end
