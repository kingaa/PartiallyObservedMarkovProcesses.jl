export dmeasure

"""
    dmeasure(object; y = obs(object), x = states(object), times=times(object), params=coef(object))

`dmeasure` is the workhorse for the evaluator of the measurement density.

The user can supply a *dmeasure* component as a function that takes data, states, parameters, and, optionally, `t`, the current time.

Calling `dmeasure` in the absence of a user-supplied *dmeasure* component results in an error.
"""
dmeasure(
    object::AbstractPompObject{T};
    y::Array{<:NamedTuple,M} = obs(object),
    x::Array{<:NamedTuple,N} = states(object),
    times::Union{T,Vector{T}} = times(object),
    params::Union{<:NamedTuple,Vector{<:NamedTuple}} = coef(object),
) where {M,N,T} = begin
    if isnothing(pomp(object).dmeasure)
        error("The *dmeasure* basic component is undefined.")
    end
    try
        times = vectorize(times)
        params = val_array(params)
        m, n, sy... = size(y)
        m1, n1, sx... = size(x)
        if m != length(times)
            error("y-times size mismatch.")
        end
        if n != length(params)
            error("y-params size mismatch.")
        end
        if (m1 != m || n1 != n)
            error("x-y size mismatch.")
        end
        if (length(sy)==0 || length(sx)==0)
            error("x and y should be arrays with at least 3 dimensions.")
        end
        x = val_array(x,m,n)
        y = val_array(y,m,n)
        if size(y,3) != 1
            error("replicate data not allowed.")
        end
        ell = [pomp(object).dmeasure(;t=times[k],y[k,j,1]...,x[k,j,i]...,params[j]...)
               for k ∈ eachindex(times), j ∈ eachindex(params), i ∈ axes(x,3)]
        reshape(ell,m,n,sx...)
    catch e
        if isa(e,UndefKeywordError)
            error("in `dmeasure`: parameter " * string(e.var) * " undefined.")
        elseif hasproperty(e,:msg)
            error("in `dmeasure`: " * e.msg)
        else
            throw(e)            # COV_EXCL_LINE
        end
    end
end
