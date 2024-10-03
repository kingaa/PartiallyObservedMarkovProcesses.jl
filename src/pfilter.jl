mutable struct PfilterdPompObject{T,X,Y} <: AbstractPompObject{T}
    pompobj::PompObject{T}
    Np::Integer
    params::NamedTuple
    eff_sample_size::Array{Float64,1}
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
        eff_sample_size = Array{Float64}(undef,length(t))
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
            cond_logLik[k],eff_sample_size[k] = pfilt_internal!(ell,p)
            xf = xp[p,:,1]
            t0 = t[k]
        end
        PfilterdPompObject{T,X,eltype(y)}(
            object,
            Np,
            params[1],
            eff_sample_size,
            cond_logLik,
            sum(cond_logLik)
        )
    catch e
        if hasproperty(e,:msg)
            error("in `pfilter`: " * e.msg)
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

pfilt_internal!(
    w::AbstractArray{<:Real,M},
    p::AbstractVector{Int64},
) where M = begin
    n = length(w)
    wmax = -Inf
    for k ∈ eachindex(w)
        if w[k] > wmax
            wmax = w[k]
        end
    end
    s = 0
    ss = 0
    for k ∈ eachindex(w)
        w[k] = w[k]/wmax
        s += w[k]
        ss += w[k]*w[k]
        w[k] = s
    end
    @assert s > 0 "sum of weights should be positive!"
    logLik = log(wmax*s/n)
    ess = s*s/ss
    du = s/length(p)
    u = -du*Random.rand()
    i = 1
    for j ∈ eachindex(p)
        u += du
        while (u > w[i] && i < n)
            i += 1
        end
        p[j] = i
    end
    logLik, ess
end
