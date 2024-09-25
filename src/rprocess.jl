export rprocess

"""
    rprocess(object; x0, t0 = timezero(object), times=times(object), params=coef(object))

`rprocess` is the workhorse for the simulator of the process

The user can supply an *rprocess* component as a function that takes states, parameters, and current time (`t`) and returns the updated time and state.

If there is no user-supplied *rprocess* component, the dynamics are trivial.
"""
rprocess(
    object::AbstractPompObject{T};
    x0::Union{X,Array{X,N}},
    t0::T = timezero(object),
    times::Union{T,Vector{T}} = times(object),
    params::Union{P,Vector{P}},
) where {T,N,X,P<:NamedTuple} = begin
    try
        times = val_array(times)
        params = val_array(params)
        x0 = val_array(x0,length(params))
        x = Array{X}(undef,size(x0)...,length(times))
        f = pomp(object).rprocess
        if isnothing(f)         # default behavior is persistence
            for k ∈ eachindex(times)
                x[:,:,k] = x0
            end
        else
            for i ∈ axes(x0,1), j ∈ eachindex(params)
                t = t0
                xx = x0[i,j]
                for k ∈ eachindex(times)
                    while t < times[k]
                        (t,xx...) = f(;t=t,xx...,params[j]...)
                    end
                    x[i,j,k] = xx
                end
            end
        end
        x
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
