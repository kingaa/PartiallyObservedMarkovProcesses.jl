mutable struct PfilterdPompObject{T,X} <: AbstractPompObject{T}
    pompobj::PompObject{T}
    Np::Integer
    params::NamedTuple
    cond_logLik::Array{Float64,1}
    logLik::Float64
end

export pomp

pomp(object::PfilterdPompObject) = object.pompobj

export coef

"""
    coef(object::PfilterdPompObject)

`coef` extracts the parameters stored in a *PfilterdPompObject*.
"""
coef(object::PfilterdPompObject) = object.params

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
        t = times(object)
        y = reshape(obs(object),1,1,length(t))
        params = val_array(params)
        xf = rinit(object,t0=t0,params=params,nsim=Np)
        X = eltype(xf)
        xp = Array{X}(undef,size(xf)...,1)
        cond_logLik = Array{Float64}(undef,length(t))
        p = Array{Int64}(undef,Np)
        ell = Array{Float64}(undef,1,Np,1,1)
        for k ∈ eachindex(t)
            rprocess!(
                object,
                xp,
                x0=xf,
                t0=t0,
                times=t[[k]],
                params=params
            )
            dmeasure!(
                object,
                ell,
                times=t[[k]],
                y=y[:,:,[k]],
                x=xp,
                params=params,
                give_log=false
            )
            cond_logLik[k] = log(mean(ell))
            systematic_resample!(p,ell)
            xf = xp[p,:,1]
            t0 = t[k]
        end
        PfilterdPompObject{T,X}(
            object,
            Np,
            params[1],
            cond_logLik,
            sum(cond_logLik)
        )
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
