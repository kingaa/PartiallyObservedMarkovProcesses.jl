export rprocess

"""
    rprocess(object; x0, t0 = timezero(object), times=times(object), params=coef(object))

`rprocess` is the workhorse for the simulator of the process

The user can supply an *rprocess* component as a function that takes states, parameters, and current time (`t`) and returns the updated time and state.

Calling `rprocess` in the absence of a user-supplied *rprocess* component results in an error.
"""
rprocess(
    object::PompObject;
    x0::Union{<:NamedTuple,Vector{<:NamedTuple},Array{<:NamedTuple,N}},
    t0::Real = timezero(object),
    times::Union{<:Real,Vector{<:Real}} = times(object),
    params::Union{<:NamedTuple,Vector{<:NamedTuple}} = coef(object),
) where N = begin
    if isnothing(object.rprocess)
        error("The *rprocess* basic component is undefined.")
    end
    try
        times = time_vector(times)
        params = val_array(params)
        x0 = val_array(x0,length(params),1)
        tx = typeof(x0)
        X = tx(undef,size(x0,1),length(params),length(times))
        t = t0
        for k ∈ eachindex(times), j = eachindex(params), i = axes(x0,1)
            while (t < times[k])
                (t,x0[i,j,1]...) = object.rprocess(;t=t,x0[i,j,1]...,params[j]...)
            end
            X[i,j,k] = x0[i,j,1]
        end
        X
    catch e
        if isa(e,UndefKeywordError)
            error("in `rprocess`: parameter " * string(e.var) * " undefined.")
        elseif hasproperty(e,:msg)
            error("in `rprocess`: " * e.msg)
        else
            throw(e)            # COV_EXCL_LINE
        end
    end
end

rprocess(object::AbstractPompObject;args...) = rprocess(pomp(object);args...)
