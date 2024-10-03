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
) where {T,P} = begin
    try
        object = pomp(object;args...)
        params = val_array(params)
        t0 = timezero(object)
        t = times(object)
        y = reshape(obs(object),1,1,length(t))
        xf = rinit(object,t0=t0,params=params,nsim=Np)
        cond_logLik = Array{Float64}(undef,length(t))
        eff_sample_size = Array{Float64}(undef,length(t))
        pfilt_internal!(
            object,
            xf,
            t0,t,y,
            params,
            eff_sample_size,
            cond_logLik
        )
        PfilterdPompObject{T,eltype(xf),eltype(y)}(
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
) where P =
    pfilter(pomp(object),Np=Np,params=params,args...)

pfilter!(
    object::PfilterdPompObject;
    params::P = coef(object),
    Np::Integer = object.Np,
    args...,
) where P = begin
    try
        object.pompobj = pomp(object;args...)
        params = val_array(params)
        t0 = timezero(object)
        t = times(object)
        y = reshape(obs(object),1,1,length(t))
        xf = rinit(object,t0=t0,params=params,nsim=Np)
        pfilt_internal!(
            object,
            xf,
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

pfilt_internal!(
    object::AbstractPompObject{T},
    xf::AbstractArray{X,2},
    t0::T,
    t::AbstractVector{T},
    y::AbstractArray{Y,3},
    params::AbstractArray{P,1},
    eff_sample_size::AbstractVector{Float64},
    cond_logLik::AbstractVector{Float64},
) where {T,X,Y,P} = begin
    Np = size(xf,1)
    xp = Array{X}(undef,size(xf)...,1)
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
        logdmeasure!(
            object,
            ell,
            times=t[[k]],
            y=y[:,:,[k]],
            x=xp,
            params=params
        )
        cond_logLik[k],eff_sample_size[k] = pfilt_step_comps!(ell,p)
        xf = xp[p,:,1]
        t0 = t[k]
    end
end

pfilt_step_comps!(
    w::AbstractArray{<:Real,M},
    p::AbstractVector{Int64},
) where M = begin
    n::Int64 = length(w)
    wmax::Float64 = -Inf
    for k ∈ eachindex(w)
        if (w[k] > wmax)
            wmax = w[k]
        end
    end
    s::Float64 = 0
    ss::Float64 = 0
    for k ∈ eachindex(w)
        v::Float64 = exp(w[k]-wmax)
        s += v
        ss += v*v
        w[k] = s
    end
    @assert s > 0 "sum of weights should be positive!"
    logLik::Float64 = wmax+log(s/n)
    ess::Float64 = s*s/ss
    du::Float64 = s/length(p)
    u::Float64 = -du*Random.rand()
    i::Int64 = 1
    for j ∈ eachindex(p)
        u += du
        while (u > w[i] && i < n)
            i += 1
        end
        p[j] = i
    end
    logLik,ess
end
