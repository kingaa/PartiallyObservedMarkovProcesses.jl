export rprocess

"""
    rprocess(object; x0, t0 = timezero(object), times=times(object), params=coef(object))

`rprocess` is the workhorse for the simulator of the process

The user can supply an *rprocess* component as a function that takes states, parameters, and current time (`t`) and returns the updated time and state.

Calling `rprocess` in the absence of a user-supplied *rprocess* component results in an error.
"""
rprocess(
    object::PompObject;
    x0::Array{<:NamedTuple,N},
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
        m, n, sx... = size(x0)
        if m != 1
            error("in `rprocess`: we must have size(x0,1)==1.")
        end
        if n != length(params)
            error("in `rprocess`: x0-params dimension mismatch.")
        end
        x0 = val_array(x0,length(params))
        tx = eltype(x0)
        X = Array{tx}(undef,length(times),length(params),size(x0,2))
        for i ∈ axes(x0,2), j ∈ eachindex(params)
            t = t0
            x = x0[j,i]
            for k ∈ eachindex(times)
                while t < times[k]
                    (t,x...) = object.rprocess(;t=t,x...,params[j]...)
                end
                X[k,j,i] = x
            end
        end
        reshape(X,length(times),length(params),sx...)
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
