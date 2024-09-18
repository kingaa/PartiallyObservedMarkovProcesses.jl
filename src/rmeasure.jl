export rmeasure

"""
    rmeasure(object; x, times=times(object), params=coef(object))

`rmeasure` is the workhorse for the simulator of the measurement distribution.

The user can supply an *rmeasure* component as a function that takes states, parameters, and, optionally, `t`, the current time.

Calling `rmeasure` in the absence of a user-supplied *rmeasure* component results in an error.
"""
rmeasure(
    object::PompObject;
    x::Union{NamedTuple,Vector{<:NamedTuple},Array{<:NamedTuple,N}},
    times::Union{<:Real,Vector{<:Real}} = times(object),
    params::Union{<:NamedTuple,Vector{<:NamedTuple}} = coef(object),
) where N = begin
    if isnothing(object.rmeasure)
        error("The *rmeasure* basic component is undefined.")
    end
    try
        times = time_vector(times)
        params = val_array(params)
        x = val_array(x,length(params),length(times))
        [object.rmeasure(;t=times[k],x[i,j,k]...,params[j]...)
         for i ∈ axes(x,1), j ∈ eachindex(params), k ∈ eachindex(times)]
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

rmeasure(object::AbstractPompObject;args...) = rmeasure(pomp(object);args...)
