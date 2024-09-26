mutable struct PfilterdPompObject{T,X} <: AbstractPompObject{T}
    pompobj::PompObject{T}
    Np::Integer
    params::NamedTuple
    x0::Array{X,2}
    states::Array{X,3}
    cond_logLik::Array{Float64,1}
    logLik::Float64
end

export pomp

pomp(object::PfilterdPompObject) = object.pompobj

export coef, states

"""
    coef(object::PfilterdPompObject)

`coef` extracts the parameters stored in a *PfilterdPompObject*.
"""
coef(object::PfilterdPompObject) = object.params

"""
    states(object::AbstractPfilterdPompObject)

`states` extracts the filtered states of a *PfilterdPompObject*.
"""
states(object::PfilterdPompObject) = object.states

export pfilter

"""
    pfilter(object; params, Np = 1, args...)

`pfilter` runs a basic particle filter.
`args...` can be used to modify or unset fields.
"""
pfilter(
    object::AbstractPompObject{T};
    params::P,
    Np::Integer = 1,
    args...,
) where {T,P<:NamedTuple} = begin
    try
        object = pomp(object;args...)
        t0 = timezero(object)
        y = obs(object)
        x0 = rinit(object,t0=t0,params=params,nsim=Np)
        t = times(object)
        X = eltype(x0)
        x = Array{X}(undef,Np,1,length(t))
        xx = val_array(x0,1,1)
        cond_logLik = Array{Float64}(undef,length(t))
        for k ∈ eachindex(t)
            xx = rprocess(object,x0=xx,t0=t0,times=t[k],params=params)
            ell = dmeasure(object,times=t[k],y=y[k],x=xx,params=params,give_log=false)
            cond_logLik[k] = log(mean(ell))
            p = systematic_resample(ell,Np)
            x[:,:,k] = xx[p,:,:]
            t0 = t[k]
        end
        PfilterdPompObject{T,X}(object,Np,params,x0,x,cond_logLik,sum(cond_logLik))
    catch e
        if hasproperty(e,:msg)              # COV_EXCL_LINE
            error("in `pfilter`: " * e.msg) # COV_EXCL_LINE
        else
            throw(e)            # COV_EXCL_LINE
        end
    end
end

pfilter(
    object::PfilterdPompObject;
    params::P = coef(object),
    Np::Integer = object.Np,
    args...,
) where {P<:NamedTuple} =
    pfilter(pomp(object),Np=Np,params=params,args...)
