export simulate, states

mutable struct SimPompObject
    pompobj::PompObject
    states::Array{<:NamedTuple,N} where N
end

convert(::Type{PompObject},object::SimPompObject) = object.pompobj

pomp(object::SimPompObject;args...) = pomp(object.pompobj;args...)
pomp!(object::SimPompObject;args...) = pomp!(object.pompobj;args...)
times(object::SimPompObject;args...) = times(object.pompobj;args...)
obs(object::SimPompObject;args...) = obs(object.pompobj;args...)
coef(object::SimPompObject;args...) = coef(object.pompobj;args...)
coef!(object::SimPompObject;args...) = coef!(object.pompobj;args...)
rinit(object::SimPompObject;args...) = rinit(object.pompobj;args...)
rmeasure(object::SimPompObject;args...) = rmeasure(object.pompobj;args...)
rprocess(object::SimPompObject;args...) = rprocess(object.pompobj;args...)
states(object::SimPompObject) = object.states

"""
    simulate(object; args...)

`simulate` simulates the POMP.
`args...` can be used to modify or unset fields.
"""
simulate(
    object::PompObject;
    args...,
) = begin
    try
        object = pomp(object;args...)
        x0 = rinit(object)
        x = rprocess(object,x0=x0)
        y = rmeasure(object,x=x)
        object.data = reshape(y,length(y))
        SimPompObject(object,x)
    catch e
        if isa(e,UndefKeywordError)
            error("in `simulate`: parameter " * string(e.var) * " undefined.")
        elseif hasproperty(e,:msg)
            error("in `simulate`: " * e.msg)
        else
            throw(e)
        end
    end
end

"""
    simulate!(object; args...)

`simulate!` simulates in place.
`args...` can be used to modify or unset fields.
"""
simulate!(
    object::SimPompObject;
    args...,
) = begin
    try
        pomp!(object;args...)
        x0 = rinit(object)
        x = rprocess(object,x0=x0)
        y = rmeasure(object,x=x)
        object.object.data = reshape(y,length(y))
        object.states = x
    catch e
        if isa(e,UndefKeywordError)
            error("in `simulate!`: parameter " * string(e.var) * " undefined.")
        elseif hasproperty(e,:msg)
            error("in `simulate!`: " * e.msg)
        else
            throw(e)
        end
    end
end
