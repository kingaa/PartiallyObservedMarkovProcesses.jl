export simulate, simulate!, states

abstract type AbstractSimPompObject <: AbstractPompObject end

mutable struct SimPompObject <: AbstractSimPompObject
    pompobj::PompObject
    states::Array{<:NamedTuple,N} where N
end

pomp(object::SimPompObject) = object.pompobj
states(object::SimPompObject) = object.states
## convert(::Type{PompObject},object::SimPompObject) = object.pompobj

"""
    simulate(object; args...)

`simulate` simulates the POMP.
`args...` can be used to modify or unset fields.
"""
simulate(
    object::AbstractPompObject;
    args...,
) = begin
    try
        object = pomp(object;args...)
        x0 = rinit(object)
        x = rprocess(object,x0=x0)
        y = rmeasure(object,x=x)
        obs!(object,reshape(y,length(y)))
        SimPompObject(pomp(object),x)
    catch e
        if hasproperty(e,:msg)
            error("in `simulate`: " * e.msg)
        else
            throw(e)            # COV_EXCL_LINE
        end
    end
end

"""
    simulate!(object; args...)

`simulate!` simulates in place.
`args...` can be used to modify or unset fields.
"""
simulate!(
    object::AbstractSimPompObject;
    args...,
) = begin
    try
        pomp!(object;args...)
        x0 = rinit(object)
        x = rprocess(object,x0=x0)
        y = rmeasure(object,x=x)
        obs!(object,reshape(y,length(y)))
        object.states=x
        object
    catch e
        if hasproperty(e,:msg)
            error("in `simulate!`: " * e.msg)
        else
            throw(e)            # COV_EXCL_LINE
        end
    end
end
