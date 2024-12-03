struct PfilterdPompObject{
    T <: Time,
    P <: NamedTuple,
    A <: Union{<:NamedTuple,Nothing},
    X0 <: Union{<:NamedTuple,Nothing},
    X <: Union{Vector{<:NamedTuple},Nothing},
    Y <: Vector{<:NamedTuple},
    } <: AbstractPompObject{T,P,A,X0,X,Y}
    pompobj::PompObject{T,P,A,X0,X,Y}
    Np::Integer
    x0::Array{<:NamedTuple,1}
    filt::Array{<:NamedTuple,2}
    pred::Array{<:NamedTuple,2}
    weights::Array{LogLik,2}
    eff_sample_size::Array{LogLik,1}
    cond_logLik::Array{LogLik,1}
    logLik::LogLik
end

export pomp

pomp(object::PfilterdPompObject) = object.pompobj

export pfilter

"""
    pfilter(object; Np = 1, params, rinit, rprocess, logmeasure, args...)

`pfilter` runs a basic particle filter.
At least the `rinit`, `rprocess`, and `logdmeasure` basic components are needed.
`args...` can be used to modify or unset additional fields.
"""
pfilter(
    object::ValidPompData;
    Np::Integer = 1,
    params::P = coef(object),
    rinit::Union{Function,Missing} = missing,
    rprocess::Union{Function,Missing} = missing,
    logdmeasure::Union{Function,Missing} = missing,
    args...,
) where {P<:NamedTuple} = let
    try
        object = pomp(
            object;
            params=params,
            rinit=rinit,
            rprocess=rprocess,
            logdmeasure=logdmeasure,
            args...,
        )
        t0 = timezero(object)
        t = times(object)
        y = obs(object)
        x0 = POMP.rinit(object;t0=t0,nsim=Np)
        X = eltype(x0)
        xf = Array{X}(undef,Np,length(t))
        xp = Array{X}(undef,Np,length(t))
        w = Array{LogLik}(undef,Np,length(t))
        cond_logLik = Array{LogLik}(undef,length(t))
        eff_sample_size = Array{LogLik}(undef,length(t))
        pfilter_internal!(
            object,
            x0,
            reshape(xf,Np,1,length(t)),
            reshape(xp,Np,1,length(t)),
            reshape(w,1,Np,1,length(t)),
            t0,t,
            reshape(y,1,1,length(t)),
            eff_sample_size,
            cond_logLik
        )
        PfilterdPompObject(
            object,
            Np,
            vec(x0),xf,xp,w,
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

"""
    pfilter(object; Np = object.Np, args...)

Running `pfilter` on a `PfilterdPompObject` re-runs the particle filter.
One can adjust the parameters, number of particles (`Np`), or pomp model components.
"""
pfilter(
    object::PfilterdPompObject;
    Np::Integer = object.Np,
    args...,
) =
    pfilter(pomp(object;args...);Np=Np)

pfilter_internal!(
    object::AbstractPompObject,
    x0::AbstractArray{X,2},
    xf::AbstractArray{X,3},
    xp::AbstractArray{X,3},
    w::AbstractArray{LogLik,4},
    t0::T,
    t::AbstractVector{T},
    y::AbstractArray{Y,3},
    eff_sample_size::AbstractVector{LogLik},
    cond_logLik::AbstractVector{LogLik},
) where {T,X,Y} = let
    for k ∈ eachindex(t)
        @inbounds rprocess!(
            object,
            @view(xp[:,:,[k]]);
            x0=x0,
            t0=t0,
            times=@view(t[[k]])
        )
        @inbounds logdmeasure!(
            object,
            @view(w[:,:,:,[k]]);
            times=@view(t[[k]]),
            y=@view(y[:,:,[k]]),
            x=@view(xp[:,:,[k]])
        )
        @inbounds pfilt_step_comps!(
            @view(cond_logLik[k]),
            @view(eff_sample_size[k]),
            @view(w[1,:,1,k]),
            @view(xp[:,:,k]),
            @view(xf[:,:,k]),
        )
        @inbounds t0 = t[k]
        @inbounds x0 = view(xf,:,:,k)
    end
end

pfilt_step_comps!(
    logLik::AbstractArray{W,0},
    ess::AbstractArray{W,0},
    w::AbstractArray{W,1},
    xp::AbstractArray{X,2},
    xf::AbstractArray{X,2},
    n::Int64 = length(w),
) where {W<:Real,X<:NamedTuple} = let
    p = Array{Int64}(undef,n)
    wmax::W = -Inf
    s::W = 0
    ss::W = 0
    for k ∈ eachindex(w)
        @inbounds wmax = (w[k] > wmax) ? w[k] : wmax
    end
    if isfinite(wmax)
        for k ∈ eachindex(w)
            @inbounds v::W = exp(w[k]-wmax)
            s += v
            ss += v*v
            @inbounds w[k] = s
        end
        ess[] = s*s/ss
        logLik[] = wmax+log(s/n)
        du::W = s/n
        u::W = -du*rand(LogLik)
        i::Int64 = 1
        for j ∈ axes(xf,1)
            u += du
            @inbounds while (u > w[i] && i < n)
                i += 1
            end
            @inbounds p[j] = i
        end
        @inbounds @views xf[:,:] = xp[p,:]
    else
        s = 0
        ss = 0
        wmax = 0
        ess[] = 0
        logLik[] = -Inf
        @inbounds @views xf[:,:] = xp[:,:]
    end
end
