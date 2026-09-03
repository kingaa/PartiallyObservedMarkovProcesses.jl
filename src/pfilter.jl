import Statistics: mean

struct PfilterdPompObject{
    T <: Time,
    X <: NamedTuple,
    P <: PompObject{T,X},
    W <: AbstractFloat
    } <: AbstractPompObject
    pompobj::P
    Np::Int
    x0::Array{X,1}
    filt::Array{X,2}
    pred::Array{X,2}
    weights::Array{W,2}
    eff_sample_size::Array{W,1}
    cond_logLik::Array{W,1}
    logLik::W
    trigger::Float64
    target::Float64
end

pomp(object::PfilterdPompObject) = object.pompobj
logLik(object::PfilterdPompObject) = object.logLik
eff_sample_size(object::PfilterdPompObject) = object.eff_sample_size
cond_logLik(object::PfilterdPompObject) = object.cond_logLik

"""
    pfilter(object; Np = 1, params, rinit, rprocess, logmeasure,
            trigger = 1, target = 0, kwargs...)

`pfilter` runs a sequential Monte Carlo computation, also known as a
particle filter.  At least the `rinit`, `rprocess`, and `logdmeasure`
basic components are needed.  Resampling is triggered whenever the
effective sample size falls below `trigger*Np`.  The resampling is
performed so that the weights are renormalized to the power `target`,
i.e., if `target = α`, `w` is a weight on call and `W` is the
corresponding renormalized weight, then `W ∝ wᵅ`.

As in other POMP.jl functions, `kwargs...` can be used to modify or
unset additional fields in the `AbstractPompObject` `object`.
"""
pfilter(
    object::ValidPompData;
    Np::Integer = 1,
    params::P = coef(object),
    rinit::Union{Function,Nothing,Missing} = missing,
    rprocess::Union{PompPlugin,Nothing,Missing} = missing,
    logdmeasure::Union{Function,Nothing,Missing} = missing,
    trigger = 1, target = 0,
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
    xf = similar(x0,length(t),Np)
    xp = similar(x0,length(t),Np)
    xt = similar(x0,length(t))
    w = Array{LogLik}(undef,length(t),Np)
    cond_logLik = similar(w,length(t))
    eff_sample_size = similar(w,length(t))
    perm = Array{Int}(undef,length(t),Np)
    trigger = Float64(clamp(trigger,0,1))
    target = Float64(clamp(target,0,1))
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
        perm,
        trigger,
        target,
    )
    i = trace_ancestry!(xt,xf,perm)
    ## FIXME: check that the first index in the ancestry is correct
    PfilterdPompObject(
        PompObject(object,init_state=x0[i],states=xt),
        Np,vec(x0),xf,xp,w,
        eff_sample_size,
        cond_logLik,
        sum(cond_logLik),
        trigger,target
    )
end

"""
    pfilter(object; Np, trigger, target, kwargs...)

Running `pfilter` on a `PfilterdPompObject` re-runs the particle filter.
One can adjust the parameters, number of particles (`Np`), or pomp model components.
"""
pfilter(
    object::PfilterdPompObject;
    Np::Integer = object.Np,
    trigger = object.trigger,
    target = object.target,
    kwargs...,
) = pfilter(pomp(object; kwargs...); Np, trigger, target)

pfilter(_...) = error("Incorrect call to `pfilter`.")

pfilter_internal!(
    object::AbstractPompObject,
    x0::AbstractArray{X,2},
    xf::AbstractArray{X,3},
    xp::AbstractArray{X,3},
    w::AbstractArray{W,4},
    t0::T,
    t::AbstractArray{T,1},
    y::AbstractArray{Y,3},
    eff_sample_size::AbstractArray{W,1},
    cond_logLik::AbstractArray{W,1},
    perm::AbstractArray{I,2},
    trigger::Float64,
    target::Float64,
) where {W<:AbstractFloat,T<:Time,X<:NamedTuple,Y<:NamedTuple,I<:Integer} = begin
    wprop = ones(W,size(x0,2))
    work = similar(wprop)
    @inbounds for k ∈ eachindex(t)
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
            wprop, work,
            @view(perm[k,:]),
            @view(xp[k,1,:]),
            @view(xf[k,1,:]),
            trigger, target,
        )
        t0 = t[k]
        x0 = view(xf,k,:,:)
    end
    nothing
end

pfilt_step_comps!(
    logLik::AbstractArray{W,0},
    ess::AbstractArray{W,0},
    logw::AbstractArray{W,1},
    w::AbstractArray{W,1},
    work::AbstractArray{W,1},
    p::AbstractArray{I,1},
    xp::AbstractArray{X,1},
    xf::AbstractArray{X,1},
    trigger::Float64,
    target::Float64,
    n::Integer = length(logw),
) where {W<:AbstractFloat,I<:Integer,X<:NamedTuple} = begin
    logwmax = compute_ess_logLik!(ess, logLik, logw, w)
    if isfinite(logwmax) && ess[] ≤ trigger*n
        systematic_resample!(p, w, work, target)
        @inbounds xf .= xp[p]
    else
        p .= collect(eachindex(p))
        xf .= xp
    end
    nothing
end

## This function computes the effective sample size (ess) and log
## likelihood (logLik).  It applies the weights in `w` to the
## log-weights in `logw`, over-writing both. On return, `w =
## exp(logw)`. It returns the maximum of `logw`. The correctness of
## this function depends on `w` having unit mean on call. It is
## guaranteed to have unit mean on exit.
compute_ess_logLik!(
    ess::AbstractArray{W,0},
    logLik::AbstractArray{W,0},
    logw::AbstractArray{W,1},
    w::AbstractArray{W,1},
) where {W <: AbstractFloat} = begin
    logwmax::W = maximum(logw)
    @assert(
        !isnan(logwmax) && logwmax < Inf,
        "invalid NaN or +∞ log likelihood in `pfilter`"
    )
    @assert length(w)==length(logw)
    if isfinite(logwmax)
        s::W = 0
        ss::W = 0
        @inbounds for k ∈ eachindex(logw)
            logw[k] += log(w[k])-logwmax
            v::W = exp(logw[k])
            s += v
            ss += v*v
            w[k] = v
        end
        lik = s/length(w)       # unit-mean assumption is needed here
        ess[] = s*s/ss
        s = log(lik)
        logLik[] = logwmax+s
        logw .-= s
        w ./= lik
    else
        ess[] = 0
        logLik[] = W(-Inf)
        logw .= zero(W)
        w .= one(W)
    end
    logwmax
end

## This function performs resampling. The indices of the selected
## particles are returned in `p`, and the weights given in `w` are
## renormalized upon return. The vector `ucum` is working memory that
## is overwritten.
systematic_resample!(
    p::AbstractArray{I,1},
    w::AbstractArray{W,1},
    ucum::AbstractArray{W,1},
    α::Float64,  # the power to which the weights will be renormalized
) where {I,W} = begin
    @assert length(ucum)==length(w)==length(p)
    s::W = 0
    β = 1-α
    @inbounds for j ∈ eachindex(w)
        s += w[j]^β
        ucum[j] = s
    end
    n::I = length(w)
    i::I = 1
    du::W = s/n
    u::W = -du*rand(W)
    @inbounds for j ∈ eachindex(p)
        u += du
        while (u > ucum[i] && i < n)
            i += 1
        end
        p[j] = i
    end
    w .^= α
    w ./= mean(w) # Other functions rely on the weights having unit mean.
    nothing
end

trace_ancestry!(
    traj::AbstractArray{X,1},
    filt::AbstractArray{X,2},
    perm::AbstractArray{I,2},
) where {X,I} = begin
    @assert size(traj,1)==size(perm,1)
    @assert size(filt)==size(perm)
    j::I = rand(axes(perm,2))
    @inbounds for i ∈ Iterators.reverse(axes(perm,1))
        traj[i] = filt[i,j]
        j = perm[i,j]
    end
    j
end

pretty_string(object::PfilterdPompObject) = begin
    pretty_string(pomp(object)) *
        ", Np=$(object.Np)" *
        ", logLik=$(round(object.logLik,digits=2))"
end
