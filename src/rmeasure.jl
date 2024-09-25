export rmeasure

"""
    rmeasure(object; x, times=times(object), params=coef(object))

`rmeasure` is the workhorse for the simulator of the measurement distribution.

The user can supply an *rmeasure* component as a function that takes states, parameters, and, optionally, `t`, the current time.
"""
rmeasure(
    object::AbstractPompObject{T};
    x::Array{X,N},
    times::Union{T,Vector{T}} = times(object),
    params::Union{P,Vector{P}},
) where {N,T,X<:NamedTuple,P<:NamedTuple} = begin
    try
        times = val_array(times)
        params = val_array(params)
        x = val_array(x,length(params),length(times))
        f = pomp(object).rmeasure
        if isnothing(f)
            Array{NamedTuple}(undef,size(x)...)
        else
            [f(;t=times[k],x[i,j,k]...,params[j]...)
             for i ∈ axes(x,1), j ∈ eachindex(params), k ∈ eachindex(times)]
        end
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
