export rmeasure

"""
    rmeasure(object; x, times=times(object), params=coef(object))

`rmeasure` is the workhorse for the simulator of the measurement distribution.

The user can supply an *rmeasure* component as a function that takes states, parameters, and, optionally, `t`, the current time.

Calling `rmeasure` in the absence of a user-supplied *rmeasure* component results in an error.
"""
rmeasure(
    object::AbstractPompObject;
    x::Array{<:NamedTuple,N},
    times::Union{<:Real,Vector{<:Real}} = times(object),
    params::Union{<:NamedTuple,Vector{<:NamedTuple}} = coef(object),
) where N = begin
    if isnothing(pomp(object).rmeasure)
        error("The *rmeasure* basic component is undefined.")
    end
    try
        times = vectorize(times)
        params = val_array(params)
        m, n, sx... = size(x)
        if length(sx)==0
            error("in `rmeasure`: x should be an array with at least 3 dimensions.")
        end
        if m != length(times)
            error("in `rmeasure`: x-times dimension mismatch.")
        end
        if n != length(params)
            error("in `rmeasure`: x-params dimension mismatch.")
        end
        x = val_array(x,m,n)
        y = [pomp(object).rmeasure(;t=times[k],x[k,j,i]...,params[j]...)
             for k ∈ eachindex(times), j ∈ eachindex(params), i ∈ axes(x,3)]
        reshape(y,m,n,sx...)
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
