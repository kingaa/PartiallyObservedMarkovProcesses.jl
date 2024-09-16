export rmeasure

"""
    rmeasure(object; state, params=object.params, time=object.time)

`rmeasure` is the workhorse for the simulator of the measurement distribution.

The user can supply an *rmeasure* component as a function that takes states, parameters, and, optionally, `t`, the current time.

Calling `rmeasure` in the absence of a user-supplied *rmeasure* component results in an error.
"""
rmeasure = function (
    object::PompObject;
    x::Union{NamedTuple,Vector{<:NamedTuple},Array{<:NamedTuple,N}},
    params::Union{NamedTuple,Vector{<:NamedTuple}} = object.params,
    time::Union{Real,Vector{Real}} = object.time
    ) where N
    if isnothing(object.rmeasure)
        error("The *rmeasure* basic component is undefined.")
    end
    try
        time = time_vector(time)
        params = val_array(params)
        x = val_array(x,length(params),length(time))
        [object.rmeasure(;t=time[k],x[i,j,k]...,params[j]...)
         for i ∈ axes(x,1), j ∈ eachindex(params), k ∈ eachindex(time)]
    catch e
        if isa(e,UndefKeywordError)
            error("in `rmeasure`: parameter " * e.var * " undefined.")
        elseif hasproperty(e,:msg)
            error("in `rmeasure`: " * e.msg)
        else
            throw(e)
        end
    end
end
