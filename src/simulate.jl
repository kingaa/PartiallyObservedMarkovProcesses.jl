mutable struct SimPompObject{T,X,Y} <: AbstractPompObject{T}
    pompobj::PompObject{T}
    x0::Array{X,2}
    states::Array{X,3}
    obsns::Array{Y,3}
    params::Vector{<:NamedTuple}
end

export pomp

pomp(object::SimPompObject) = object.pompobj

export coef, states, obs

"""
    coef(object::SimPompObject)

`coef` extracts the parameters stored in a *SimPompObject*.
"""
coef(object::SimPompObject) = object.params

"""
    states(object::AbstractSimPompObject)

`states` extracts the state trajectory of a *SimPompObject*.
"""
states(object::SimPompObject) = object.states

"""
    obs(object::AbstractSimPompObject)

`obs` extracts the observations of a *SimPompObject*.
"""
obs(object::SimPompObject) = object.obsns

export simulate, simulate!

"""
    simulate(object; params, nsim = 1, args...)

`simulate` simulates the POMP.
`args...` can be used to modify or unset fields.
"""
simulate(
    object::AbstractPompObject{T};
    params::Union{P,Vector{P}},
    nsim::Integer = 1,
    args...,
) where {T,P<:NamedTuple} = begin
    try
        object = pomp(object;args...)
        params = val_array(params)
        x0 = rinit(object,params=params,nsim=nsim)
        x = rprocess(object,x0=x0,params=params)
        y = rmeasure(object,x=x,params=params)
        X = eltype(x)
        Y = eltype(y)
        SimPompObject{T,X,Y}(
            object,
            x0,
            x,
            y,
            params
        )
    catch e
        if hasproperty(e,:msg)               # COV_EXCL_LINE
            error("in `simulate`: " * e.msg) # COV_EXCL_LINE
        else
            throw(e)            # COV_EXCL_LINE
        end
    end
end

simulate(
    object::SimPompObject{T,X,Y};
    params::Union{P,Vector{P}} = coef(object),
    nsim::Integer = 1,
    args...,
) where {T,X,Y,P<:NamedTuple} =
    simulate(pomp(object),params=params,nsim=nsim,args...)

"""
    simulate!(object; args...)

`simulate!` simulates in place.
`args...` can be used to modify or unset fields.
"""
simulate!(
    object::SimPompObject;
    params::Union{P,Vector{P}} = coef(object),
    nsim::Integer = 1,
    args...,
) where {P<:NamedTuple} = begin
    try
        pomp!(pomp(object);args...)
        x0 = rinit(object,params=params,nsim=nsim)
        x = rprocess(object,x0=x0,params=params)
        y = rmeasure(object,x=x,params=params)
        object.x0 = x0
        object.states = x
        object.obsns = y
        object.params = params
        object
    catch e
        if hasproperty(e,:msg)               # COV_EXCL_LINE
            error("in `simulate`: " * e.msg) # COV_EXCL_LINE
        else
            throw(e)            # COV_EXCL_LINE
        end
    end
end
