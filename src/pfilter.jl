import Statistics: mean

"""
    PfilterdPompObject

Created by a call to [`pfilter`](@ref `pfilter`), this struct holds
the results of the particle-filter computation.  In particular, one
can re-run a `pfilter` computation, optionally with modifications to
the model, the model parameters, or algorithm settings, by calling
`pfilter` on a `PfilterdPompObject`.
"""
struct PfilterdPompObject{
    T <: Time,
    X <: NamedTuple,
    Y <: NamedTuple,
    P <: PompObject{T,X,Y},
    W <: AbstractFloat
    } <: AbstractPompObject
    "pomp object"
    pompobj::P
    "number of particles"
    Np::Int
    "initial states"
    x0::Array{X,1}
    "filter distribution"
    filt::Array{X,2}
    "prediction distribution"
    pred::Array{X,2}
    "Boolean indicator of whether resampling has been performed"
    resample::Array{Bool,1}
    "log particle weights (before resampling)"
    logweights::Array{W,2}
    "final weights"
    weights::Array{W,1}
    "effective sample size"
    eff_sample_size::Array{W,1}
    "conditional log likelihoods"
    cond_logLik::Array{W,1}
    "sample-size fraction that triggers resampling; missing value is equivalent to 1.0"
    trigger::Union{Float64,Missing}
    "renormalization power; missing value is equivalent to 0.0"
    target::Union{Float64,Missing}
    "log likelihood estimate (=sum of `cond_logLik`)"
    logLik::W
end

pomp(object::PfilterdPompObject) = object.pompobj
logLik(object::PfilterdPompObject) = object.logLik
eff_sample_size(object::PfilterdPompObject) = object.eff_sample_size
cond_logLik(object::PfilterdPompObject) = object.cond_logLik

"""
    pfilter(object; Np = 1, params, rinit, rprocess, logmeasure,
            trigger, target, kwargs...)

`pfilter` runs a sequential Monte Carlo computation, also known as a
particle filter.  At least the `rinit`, `rprocess`, and `logdmeasure`
basic components are needed.  Resampling is triggered whenever the
effective sample size falls below `trigger*Np`.  The resampling is
performed so that the weights are renormalized to the power `target`,
i.e., if `target = β`, `w` is a particle weight, and `W` is the
corresponding renormalized weight, then `W ∝ wᵝ`.  One must have 0 ≤
`trigger` ≤ 1 and 0 ≤ `target` < 1.  `trigger = missing` is a synonym
for `trigger = 1` and `target = missing` is equivalent to `target =
0`.

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
    trigger::Union{Real,Missing} = missing,
    target::Union{Real,Missing} = missing,
    kwargs...,
) where {P<:NamedTuple} = begin
    object = pomp(
        object;
        params, rinit, rprocess, logdmeasure,
        kwargs...,
    )
    trigger, target = proc_trig_targ(trigger, target)
    x0 = POMP.rinit(object; nsim=Np)
    xf, xp, logw, w, cll, ess, perm, resamp = pfilter_internal!(
        object, x0,
        trigger, target,
    )
    xt = similar(x0, length(times(object)))
    i = trace_ancestry!(xt, xf, perm, w)
    PfilterdPompObject(
        PompObject(object, init_state=x0[i], states=xt),
        Np, vec(x0) ,xf, xp, resamp, logw, w,
        ess, cll,
        trigger, target,
        sum(cll),
    )
end

"""
    pfilter(object; Np, trigger, target, kwargs...)

Running `pfilter` on a `PfilterdPompObject` re-runs the particle
filter.  Additional arguments adjust the number of particles (`Np`),
`trigger` and/or `target` specifications, model parameters, or basic
model components.
"""
pfilter(
    object::PfilterdPompObject;
    Np::Integer = object.Np,
    trigger::Union{Real,Missing} = object.trigger,
    target::Union{Real,Missing} = object.target,
    kwargs...,
) = pfilter(pomp(object); Np, trigger, target, kwargs...)

pfilter(_...) = error("Incorrect call to `pfilter`.")

## allocate memory and run main loop
pfilter_internal!(
    object, x0, args...,
) = begin
    t0 = timezero(object)
    t = times(object)
    y = obs(object)
    Np = length(x0)
    xf = similar(x0,length(t),Np) # filter distribution
    xp = similar(x0,length(t),Np) # prediction distribution
    logw = Array{LogLik}(undef,length(t),Np) # log weights
    cll = similar(logw,length(t)) # conditional log likelihood
    ess = similar(logw,length(t)) # effective sample size
    perm = Array{Int}(undef,length(t),Np) # sampled indices
    resamp = Array{Bool}(undef,length(t)) # indicator of resampling
    w = ones(LogLik,Np)
    pfilter_loop!(
        object,
        t0, t, x0,
        reshape(xf,length(t),1,Np),
        reshape(xp,length(t),1,Np),
        reshape(y,length(t),1,1),
        reshape(logw,length(t),1,Np,1),
        w, ess, cll, perm, resamp,
        args...,
    )
    xf, xp, logw, w, cll, ess, perm, resamp
end

## main loop for weighted particle filter
pfilter_loop!(
    object::AbstractPompObject,
    t0::T,
    t::AbstractArray{T,1},
    x0::AbstractArray{X,2},
    xf::AbstractArray{X,3},
    xp::AbstractArray{X,3},
    y::AbstractArray{Y,3},
    logw::AbstractArray{W,4},
    wprop::AbstractArray{W,1},
    eff_sample_size::AbstractArray{W,1},
    cond_logLik::AbstractArray{W,1},
    perm::AbstractArray{I,2},
    resample::AbstractArray{Bool,1},
    trigger::Float64,
    target::Float64,
) where {T<:Time,X<:NamedTuple,W<:AbstractFloat,Y<:NamedTuple,I<:Integer} = begin
    work = similar(wprop)
    for k ∈ eachindex(t)
        pfilter_step!(
            object, k, t0, t, x0, xp, xf, y,
            logw, cond_logLik, eff_sample_size,
            perm, resample, work,
            wprop, trigger, target,
        )
        t0 = t[k]
        x0 = view(xf,k,:,:)
    end
    nothing
end

## main loop for unweighted particle filter
pfilter_loop!(
    object::AbstractPompObject,
    t0::T,
    t::AbstractArray{T,1},
    x0::AbstractArray{X,2},
    xf::AbstractArray{X,3},
    xp::AbstractArray{X,3},
    y::AbstractArray{Y,3},
    logw::AbstractArray{W,4},
    _::AbstractArray{W,1},
    eff_sample_size::AbstractArray{W,1},
    cond_logLik::AbstractArray{W,1},
    perm::AbstractArray{I,2},
    resample::AbstractArray{Bool,1},
    trigger::Missing,
    target::Missing,
) where {T<:Time,X<:NamedTuple,W<:AbstractFloat,Y<:NamedTuple,I<:Integer} = begin
    work = Array{W}(undef,size(x0,2))
    for k ∈ eachindex(t)
        pfilter_step!(
            object, k, t0, t, x0, xp, xf, y,
            logw, cond_logLik, eff_sample_size,
            perm, resample, work,
        )
        t0 = t[k]
        x0 = view(xf,k,:,:)
    end
    nothing
end

## one step in particle-filter loop
pfilter_step!(
    object::AbstractPompObject,
    k::Integer,
    t0::T,
    t::AbstractArray{T,1},
    x0::AbstractArray{X,2},
    xp::AbstractArray{X,3},
    xf::AbstractArray{X,3},
    y::AbstractArray{Y,3},
    logw::AbstractArray{W,4},
    cond_logLik::AbstractArray{W,1},
    eff_sample_size::AbstractArray{W,1},
    perm::AbstractArray{I,2},
    resample::AbstractArray{Bool,1},
    args...,
) where {W<:AbstractFloat,T<:Time,X<:NamedTuple,Y<:NamedTuple,I<:Integer} = begin
    advance_particles!(
        object, t0, @view(t[[k]]),
        x0, @view(xp[[k],:,:]),
        @view(y[[k],:,:]),
        @view(logw[[k],:,:,:]),
    )
    pfilt_step_comps!(
        @view(cond_logLik[k]),
        @view(eff_sample_size[k]),
        @view(logw[k,1,:,1]),
        @view(perm[k,:]),
        @view(xp[k,1,:]),
        @view(xf[k,1,:]),
        @view(resample[k]),
        args...,
    )
    nothing
end

advance_particles!(object, t0, times, x0, x, y, w) = begin
    rprocess!(object,x;x0,t0,times)
    logdmeasure!(object,w;times,y,x)
    nothing
end

pfilt_step_comps!(
    logLik::AbstractArray{W,0},
    ess::AbstractArray{W,0},
    logw::AbstractArray{W,1},
    p::AbstractArray{I,1},
    xp::AbstractArray{X,1},
    xf::AbstractArray{X,1},
    resample::AbstractArray{Bool,0},
    work::AbstractArray{W,1},
    w::AbstractArray{W,1},
    trigger::Float64,
    target::Float64,
    n::Integer = length(logw),
) where {W<:AbstractFloat,I<:Integer,X<:NamedTuple} = begin
    logwmax = compute_ess_logLik!(ess, logLik, logw, w)
    if isfinite(logwmax) && ess[] ≤ trigger*n
        systematic_resample!(p, w, work, target)
        resample[] = true
        xf .= xp[p]
    else
        p .= collect(eachindex(p))
        resample[] = false
        xf .= xp
    end
    nothing
end

pfilt_step_comps!(
    logLik::AbstractArray{W,0},
    ess::AbstractArray{W,0},
    logw::AbstractArray{W,1},
    p::AbstractArray{I,1},
    xp::AbstractArray{X,1},
    xf::AbstractArray{X,1},
    resample::AbstractArray{Bool,0},
    work::AbstractArray{W,1},
) where {W<:AbstractFloat,I<:Integer,X<:NamedTuple} = begin
    logwmax = compute_ess_logLik!(ess, logLik, logw)
    if isfinite(logwmax)
        systematic_resample!(p, logw, work)
        resample[] = true
        xf .= xp[p]
    else
        p .= collect(eachindex(p))
        resample[] = false
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
        for k ∈ eachindex(logw)
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
        w ./= lik               # enforces unit-mean on return
    else
        ess[] = 0
        logLik[] = W(-Inf)
        logw .= zero(W)
        w .= one(W)
    end
    logwmax
end

compute_ess_logLik!(
    ess::AbstractArray{W,0},
    logLik::AbstractArray{W,0},
    logw::AbstractArray{W,1},
) where {W <: AbstractFloat} = begin
    logwmax::W = maximum(logw)
    @assert(
        !isnan(logwmax) && logwmax < Inf,
        "invalid NaN or +∞ log likelihood in `pfilter`"
    )
    if isfinite(logwmax)
        s::W = 0
        ss::W = 0
        for k ∈ eachindex(logw)
            logw[k] -= logwmax
            v::W = exp(logw[k])
            s += v
            ss += v*v
        end
        lik = s/length(logw)
        ess[] = s*s/ss
        s = log(lik)
        logLik[] = logwmax+s
        logw .-= s
    else
        ess[] = 0
        logLik[] = W(-Inf)
        logw .= zero(W)
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
    β::Float64,  # the power to which the weights will be renormalized
) where {I,W} = begin
    @assert length(ucum)==length(w)==length(p)
    s::W = 0
    α = 1-β
    for j ∈ eachindex(w)
        s += w[j]^α
        ucum[j] = s
    end
    n::I = length(w)
    i::I = 1
    du::W = s/n
    u::W = -du*rand(W)
    for j ∈ eachindex(p)
        u += du
        while (u > ucum[i] && i < n)
            i += 1
        end
        p[j] = i
    end
    n = 0
    for j ∈ eachindex(p)
        if n ≠ p[j]
            n = p[j]
            s = w[n]^β
        end
        ucum[j] = s
    end
    w .= ucum
    w ./= mean(w) # Other functions rely on the weights having unit mean.
    nothing
end

systematic_resample!(
    p::AbstractArray{I,1},
    logw::AbstractArray{W,1},
    ucum::AbstractArray{W,1},
) where {I,W} = begin
    @assert length(ucum)==length(logw)==length(p)
    s::W = 0
    for j ∈ eachindex(logw)
        s += exp(logw[j])
        ucum[j] = s
    end
    n::I = length(logw)
    i::I = 1
    du::W = s/n
    u::W = -du*rand(W)
    for j ∈ eachindex(p)
        u += du
        while (u > ucum[i] && i < n)
            i += 1
        end
        p[j] = i
    end
    nothing
end

trace_ancestry!(
    traj::AbstractArray{X,1},
    filt::AbstractArray{X,2},
    perm::AbstractArray{I,2},
    weights::AbstractArray{W,1},
) where {X,I,W} = begin
    @assert size(traj,1)==size(perm,1)
    @assert size(weights,1)==size(perm,2)
    @assert size(filt)==size(perm)
    r::W = length(weights)*rand(W) ## this relies on mean(weights)=1
    j::I = 1                       ## choose a random particle
    while r > weights[j] && j < length(weights)
        r -= weights[j]
        j += 1
    end
    for i ∈ Iterators.reverse(axes(perm,1))
        traj[i] = filt[i,j]
        j = perm[i,j]
    end
    j
end

## process the trigger and target arguments
proc_trig_targ(trigger, target) = begin
    if ismissing(trigger) && !ismissing(target)
        trigger = one(Float64)
        target = Float64(target)
    elseif !ismissing(trigger) && ismissing(target)
        trigger = Float64(trigger)
        target = zero(Float64)
    elseif !ismissing(trigger) && !ismissing(target)
        trigger = Float64(trigger)
        target = Float64(target)
    end
    @assert ismissing(trigger) || 0.0 ≤ trigger ≤ 1.0 "`trigger` should be in [0,1] or missing."
    @assert ismissing(target) || 0.0 ≤ target < 1.0 "`target` should be in [0,1) or missing."
    trigger, target
end

pretty_string(object::PfilterdPompObject) = begin
    pretty_string(pomp(object)) *
        ", Np=$(object.Np)" *
        ", logLik=$(round(object.logLik,digits=2))"
end
