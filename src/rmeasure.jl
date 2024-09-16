export rmeasure

"""
    rmeasure(object; state, params=object.params, time=object.time)

`rmeasure` is the workhorse for the simulator of the measurement distribution.

The user can supply an *rmeasure* component as a function that takes states, parameters, and, optionally, `t`, the current time.

Calling `rmeasure` in the absence of a user-supplied *rmeasure* component results in an error.
"""
rmeasure = function (
    object::PompObject;
    state::Union{NamedTuple,Array{<:NamedTuple,3}},
    params::Union{NamedTuple,Vector{<:NamedTuple}} = object.params,
    time::Union{Real,Vector{Real}} = object.time
    )
    if isnothing(object.rmeasure)
        error("The *rmeasure* basic component is undefined.")
    end
    try
        time = time_vector(time)
        params = param_vector(params)
        state = state_array(state,length(params),length(time))
        [object.rmeasure(;t=time[k],state[i,j,k]...,params[j]...)
         for i ∈ axes(state,1), j ∈ eachindex(params), k ∈ eachindex(time)]
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
