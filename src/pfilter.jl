mutable struct PfilterdPompObject{T,X,Y} <: AbstractPompObject{T}
    pompobj::PompObject{T}
    Np::Integer
    params::NamedTuple
    x0::Array{X,2}
    filt::Array{X,3}
    pred::Array{X,3}
    weights::Array{Float64,3}
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

export pfilter, pfilter!

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
) where {T,P} = let
    try
        object = pomp(object;args...)
        params = val_array(params)
        t0 = timezero(object)
        t = times(object)
        y = reshape(obs(object),1,1,length(t))
        x0 = rinit(object,t0=t0,params=params,nsim=Np)
        X = eltype(x0)
        xf = Array{X}(undef,Np,1,length(t))
        xp = Array{X}(undef,Np,1,length(t))
        w = Array{Float64}(undef,Np,1,length(t))
        cond_logLik = Array{Float64}(undef,length(t))
        eff_sample_size = Array{Float64}(undef,length(t))
        pfilter_internal!(
            object,
            x0,xf,xp,w,
            t0,t,y,
            params,
            eff_sample_size,
            cond_logLik
        )
        PfilterdPompObject{T,X,eltype(y)}(
            object,
            Np,
            params[1],
            x0,xf,xp,w,
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
) where P =
    pfilter(pomp(object),Np=Np,params=params,args...)

pfilter!(
    object::PfilterdPompObject{T,X,Y};
    params::P = coef(object),
    _...,
) where {P,T,X,Y} = begin
    try
        params = val_array(params)
        t0 = timezero(object)
        t = times(object)
        y = reshape(obs(object),1,1,length(t))
        rinit!(object,object.x0,t0=t0,params=params)
        pfilter_internal!(
            object,
            object.x0,
            object.filt,
            object.pred,
            object.weights,
            t0,t,y,
            params,
            object.eff_sample_size,
            object.cond_logLik
        )
        object.logLik = sum(object.cond_logLik)
        object
    catch e
        if hasproperty(e,:msg)
            error("in `pfilter!`: " * e.msg)
        else
            throw(e)            # COV_EXCL_LINE
        end
    end
end

pfilter_internal!(
    object::AbstractPompObject{T},
    x0::AbstractArray{X,2},
    xf::AbstractArray{X,3},
    xp::AbstractArray{X,3},
    w::AbstractArray{Float64,3},
    t0::T,
    t::AbstractVector{T},
    y::AbstractArray{Y,3},
    params::AbstractArray{P,1},
    eff_sample_size::AbstractVector{Float64},
    cond_logLik::AbstractVector{Float64},
) where {T,X,Y,P} = let
    Np = size(x0,1)
    w = reshape(w,1,Np,1,length(t))
    for k ∈ eachindex(t)
        rprocess!(
            object,
            @view(xp[:,:,[k]]),
            x0=x0,
            t0=t0,
            times=@view(t[[k]]),
            params=params
        )
        logdmeasure!(
            object,
            @view(w[:,:,:,[k]]),
            times=@view(t[[k]]),
            y=@view(y[:,:,[k]]),
            x=@view(xp[:,:,[k]]),
            params=params
        )
        pfilt_step_comps!(
            @view(cond_logLik[k]),
            @view(eff_sample_size[k]),
            @view(w[1,:,1,k]),
            @view(xp[:,:,k]),
            @view(xf[:,:,k]),
        )
        t0 = t[k]
        x0 = view(xf,:,:,k)
    end
end

pfilt_step_comps!(
    logLik::AbstractArray{W,0},
    ess::AbstractArray{W,0},
    w::AbstractArray{W,1},
    xp::AbstractArray{X,2},
    xf::AbstractArray{X,2},
) where {W<:Real,X<:NamedTuple} = let
    n::Int64 = length(w)
    p = Array{Int64}(undef,n)
    wmax::W = -Inf
    for k ∈ eachindex(w)
        if (w[k] > wmax)
            wmax = w[k]
        end
    end
    s::W = 0
    ss::W = 0
    if isfinite(wmax)
        for k ∈ eachindex(w)
            v::W = exp(w[k]-wmax)
            s += v
            ss += v*v
            w[k] = s
        end
        ess[] = s*s/ss
        logLik[] = wmax+log(s/n)
        du::W = s/n
        u::W = -du*Random.rand()
        i::Int64 = 1
        for j ∈ axes(xf,1)
            u += du
            while (u > w[i] && i < n)
                i += 1
            end
            p[j] = i
        end
        @views xf[:,:] = xp[p,:]
    else
        s = 0
        ss = 0
        wmax = 0
        ess[] = 0
        logLik[] = -Inf
        @views xf[:,:] = xp[:,:]
    end
end
