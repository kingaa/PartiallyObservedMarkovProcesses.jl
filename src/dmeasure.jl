export dmeasure

"""
    dmeasure(object; times=times(object), y = obs(object), x = states(object), params)

`dmeasure` is the workhorse for the evaluator of the log measurement density.

The user can supply a *dmeasure* component as a function that takes data, states, parameters, and, optionally, `t`, the current time.
"""
dmeasure(
    object::AbstractPompObject{T};
    times::Union{T,Vector{T}} = times(object),
    y::Array{Y,M} = obs(object),
    x::Array{X,N},
    params::Union{P,Vector{P}},
) where {M,N,T,Y<:NamedTuple,X<:NamedTuple,P<:NamedTuple} = begin
    try
        times = val_array(times)
        params = val_array(params)
        if length(x)!=length(y)
            error("x and y should be of the same size.")
        end
        x = val_array(x,length(params),length(times))
        y = val_array(y,length(params),length(times))
        f = pomp(object).dmeasure
        if isnothing(f)         # default behavior is no information
            zeros(Float64,size(x)...)
        else
            [f(;t=times[k],y[i,j,k]...,x[i,j,k]...,params[j]...)
             for i ∈ axes(x,1), j ∈ eachindex(params), k ∈ eachindex(times)]
        end
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
