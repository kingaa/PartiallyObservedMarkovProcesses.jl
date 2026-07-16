struct PfilterdPompObject{
    P <:PompObject,
    Z <: NamedTuple,
    W <: AbstractFloat
    } <: AbstractPompObject
    pompobj::P
    Np::Int
    x0::Array{Z,1}
    filt::Array{Z,2}
    pred::Array{Z,2}
    traj::Array{Z,1}
    weights::Array{W,2}
    eff_sample_size::Array{W,1}
    cond_logLik::Array{W,1}
    logLik::W
end


pomp(object::PfilterdPompObject) = object.pompobj

logLik(object::PfilterdPompObject) = object.logLik
eff_sample_size(object::PfilterdPompObject) = object.eff_sample_size
cond_logLik(object::PfilterdPompObject) = object.cond_logLik


"""
    pfilter(object; Np = 1, params, rinit, rprocess, logmeasure, kwargs...)

`pfilter` runs a basic particle filter.
At least the `rinit`, `rprocess`, and `logdmeasure` basic components are needed.
`kwargs...` can be used to modify or unset additional fields.
"""
pfilter(
    object::ValidPompData;
    Np::Integer = 1,
    params::P = coef(object),
    rinit::Union{Function,Nothing,Missing} = missing,
    rprocess::Union{PompPlugin,Nothing,Missing} = missing,
    logdmeasure::Union{Function,Nothing,Missing} = missing,
    kwargs...,
) where {P<:NamedTuple} = begin
    object = pomp(
        object;
        params,rinit,rprocess,logdmeasure,
        kwargs...,
    )
    t0 = timezero(object)
    t = times(object)
    y = obs(object)
    x0 = POMP.rinit(object;t0,nsim=Np)
    @assert(size(x0)==(1,Np))
    xf = similar(x0,length(t),Np)
    xp = similar(x0,length(t),Np)
    xt = similar(x0,length(t))
    w = similar(Array{LogLik},length(t),Np)
    cond_logLik = similar(w,length(t))
    eff_sample_size = similar(w,length(t))
    perm = similar(Array{Int64},length(t),Np)
    pfilter_internal!(
        object,
        x0,
        reshape(xf,length(t),1,Np),
        reshape(xp,length(t),1,Np),
        reshape(w,length(t),1,Np,1),
        t0,t,
        reshape(y,length(t),1,1),
        eff_sample_size,
        cond_logLik,
        perm
    )
    trace_ancestry!(xt,xf,perm)
    PfilterdPompObject(
        object,Np,
        vec(x0),xf,xp,xt,w,
        eff_sample_size,
        cond_logLik,
        sum(cond_logLik)
    )
end

"""
    pfilter(object; Np = object.Np, kwargs...)

Running `pfilter` on a `PfilterdPompObject` re-runs the particle filter.
One can adjust the parameters, number of particles (`Np`), or pomp model components.
"""
pfilter(
    object::PfilterdPompObject;
    Np::Integer = object.Np,
    kwargs...,
) = pfilter(pomp(object; kwargs...); Np)

pfilter(_...) = error("Incorrect call to `pfilter`.")

pfilter_internal!(
    object::AbstractPompObject,
    x0::AbstractArray{X,2},
    xf::AbstractArray{X,3},
    xp::AbstractArray{X,3},
    w::AbstractArray{LogLik,4},
    t0::T,
    t::AbstractArray{T,1},
    y::AbstractArray{Y,3},
    eff_sample_size::AbstractArray{W,1},
    cond_logLik::AbstractArray{W,1},
    perm::AbstractArray{Int,2}
) where {W<:AbstractFloat,T<:Time,X<:NamedTuple,Y<:NamedTuple} = begin
    for k ∈ eachindex(t)
        rprocess!(
            object,
            @view(xp[[k],:,:]);
            x0=x0,
            t0=t0,
            times=@view(t[[k]])
        )
        logdmeasure!(
            object,
            @view(w[[k],:,:,:]);
            times=@view(t[[k]]),
            y=@view(y[[k],:,:]),
            x=@view(xp[[k],:,:])
        )
        pfilt_step_comps!(
            @view(cond_logLik[k]),
            @view(eff_sample_size[k]),
            @view(w[k,1,:,1]),
            @view(perm[k,:]),
            @view(xp[k,1,:]),
            @view(xf[k,1,:]),
        )
        t0 = t[k]
        x0 = view(xf,k,:,:)
    end
    nothing
end

pfilt_step_comps!(
    logLik::AbstractArray{W,0},
    ess::AbstractArray{W,0},
    w::AbstractArray{W,1},
    p::AbstractArray{I,1},
    xp::AbstractArray{X,1},
    xf::AbstractArray{X,1},
    n::Integer = length(w),
) where {W<:AbstractFloat,I<:Integer,X<:NamedTuple} = begin
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
        i::I = 1
        for j ∈ eachindex(p)
            u += du
            @inbounds while (u > w[i] && i < n)
                i += 1
            end
            @inbounds p[j] = i
        end
        @inbounds xf[:] = xp[p]
    else
        s = 0
        ss = 0
        wmax = 0
        ess[] = 0
        logLik[] = W(-Inf)
        p[:] = collect(eachindex(p))
        @inbounds xf[:] = xp[:]
    end
    nothing
end

trace_ancestry!(
    traj::AbstractArray{X,1},
    filt::AbstractArray{X,2},
    perm::AbstractArray{I,2},
) where {X,I<:Integer} = begin
    @assert size(traj,1)==size(perm,1)
    @assert size(filt)==size(perm)
    j = rand(axes(perm,2))
    for i ∈ Iterators.reverse(axes(perm,1))
        @inbounds traj[i] = filt[i,j]
        @inbounds j = perm[i,j]
    end
    nothing
end

pretty_string(object::PfilterdPompObject) = begin
    pretty_string(pomp(object)) *
        ", Np=$(object.Np)" *
        ", logLik=$(round(object.logLik,digits=2))"
end
