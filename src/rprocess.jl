export rprocess

"""
    rprocess(object; params=object.params, state, t0 = object.t0, time=object.time)

`rprocess` is the workhorse for the simulator of the process

The user can supply an *rprocess* component as a function that takes states, parameters, and two times t1, t2.

Calling `rprocess` in the absence of a user-supplied *rprocess* component results in an error.
"""
rprocess = function (
    object::PompObject;
    t0::Real = object.t0,
    time::Union{Real,Vector{Real}} = object.time,
    x0::Union{NamedTuple,Vector{<:NamedTuple},Array{<:NamedTuple,N}},
    params::Union{NamedTuple,Vector{<:NamedTuple}} = object.params
    ) where N
    if isnothing(object.rprocess)
        error("The *rprocess* basic component is undefined.")
    end
    try
        time = time_vector(time)
        params = val_array(params)
        x0 = val_array(x0,length(params),1)
        tx = typeof(x0)
        X = tx(undef,size(x0,1),length(params),length(time))
        t = t0
        for k ∈ eachindex(time), j = eachindex(params), i = axes(x0,1)
            while (t < time[k])
                (t,x0[i,j,1]...) = object.rprocess(;t=t,x0[i,j,1]...,params[j]...)
            end
            X[i,j,k] = x0[i,j,1]
        end
        X
    catch e
        if isa(e,UndefKeywordError)
            error("in `rprocess`: parameter " * e.var * " undefined.")
        elseif hasproperty(e,:msg)
            error("in `rprocess`: " * e.msg)
        else
            throw(e)
        end
    end
end
