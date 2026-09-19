import Statistics: mean
import DataFrames: insertcols!

"""
    MifdPompObject

Created by a call to [`mif`](@ref `mif`), this struct holds the
results of the iterated filtering computation.  In particular, one can
re-run a `mif` computation, optionally with modifications to the
model, the model parameters, or algorithm settings, by calling `mif`
on a `MifdPompObject`.
"""
struct MifdPompObject{
    T <: Time,
    X <: NamedTuple,
    Y <: NamedTuple,
    P <: NamedTuple,
    W <: AbstractFloat,
    Q <: PfilterdPompObject{T,X,Y},
    } <: AbstractPompObject
    "PfilterdPompObject"
    pfobj::Q
    "number of iterations"
    Nmif::Int
    "cooling schedule function"
    cooling::Function
    "perturbations function"
    perturbations::Function
    "traces"
    trace::Array{P,1}
    "logLik trace"
    logLik::Array{W,1}
end

pomp(object::MifdPompObject) = pomp(object.pfobj)
eff_sample_size(object::MifdPompObject) = eff_sample_size(object.pfobj)
cond_logLik(object::MifdPompObject) = cond_logLik(object.pfobj)

"""
    logLik(object::MifdPompObject)

Returns the esimated log likelihood obtained using a [`pfilter`](@ref `pfilter`) computation after the final `mif` iteration.
"""
logLik(object::MifdPompObject) = logLik(object.pfobj)

"""
    traces(object)

Returns a `DataFrame` containing the traces of an iterated filtering
computation. Specifically, these are the iteration-by-iteration
trajectory of the model parameters and the estimated mif likelihood.
"""
traces(object::MifdPompObject) = begin
    insertcols!(melt(object.trace, :iteration), :logLik => [object.logLik...,missing])
end

"""
    pfilter(object; kwargs...)

Calling `pfilter` on the result of a [`mif`](@ref `mif`) computation
runs a particle filter.  By default, the parameters `Np`, `trigger`,
and `target` used in the `mif` computation are re-used; one can
optionally modify these.
"""
pfilter(object::MifdPompObject; kwargs...,) = pfilter(object.pfobj; kwargs...)

"""
    mif(object; Np = 1, Nmif = 1, perturbations, cooling, trigger,
        target, params, rinit, rprocess, logdmeasure, kwargs...)

Iterated filtering.  In addition to the components needed for a
[`pfilter`](@ref `pfilter`) (i.e., `Np`, `trigger`, `target`), one
must specify a perturbations function, cooling schedule, and number of
iterations.  After performing the requested `Nmif` iterations, `mif`
runs a particle filter using the estimated parameters.

## Arguments

- `object`: A `DataFrame`, `PompObject`, or vector of data.
- `Np`: number of particles to use in the filtering.
- `Nmif`: number of MIF iterations to perform.
- `perturbations`: a function that returns perturbed versions of some
  or all of the model parameters. See below for details.
- `cooling`: a function that specifies the MIF cooling schedule. See
  below for details.
- `trigger`, `target`: see [`pfilter`](@ref `pfilter`).
- `params`: `NamedTuple` of model parameters.
- `rinit`, `rprocess`, `logdmeasure`: necessary basic model components.
- `kwargs...`: other arguments are passed to [`pomp`](@ref `pomp`).

## Return value

`mif` returns a `MifdPompObject`.  This contains the
`PfilterdPompObject` containing the results of the final particle
filter computation.  It also records the algorithmic parameters (i.e.,
`Np`, `Nmif`, `perturbations` and `cooling` functions,
`trigger` and `target`).

## Perturbations

Iterated filtering works by applying a random perturbation to some or
all of the model parameters during a particle-filter computation. The
`perturbations` argument specifies the nature of the perturbations to
be applied.  Specifically, it should be function with signature
`f(scale, lag; params...)` where `scale` is the relative scale of the
perturbations (a fraction between 0 and 1) and `lag` is an integer
indicating the observation number (ranging from 1 to `n` if there are
`n` observations); at the initial time, `lag = 0`.  The named
arguments `params` should be the model parameters to be perturbed.  In
particular, when the function is called, these arguments will contain
the values of the corresponding model parameters.

The function should return a `NamedTuple` containing the perturbed
parameters.  Thus, for example, if one is attempting to estimate
parameters `α` and `β`, while leaving parameters `γ` and `δ` fixed,
one might furnish a function such as the following as the
`perturbations` argument to `mif`:
```
    p(scale, lag; α, β, _...) = begin
        α=rand(LogNormal(log(α),scale*0.02))
        β=rand(LogNormal(log(β),scale*0.02))
        (;α,β)
    end
```
Note that this function allows for, but ignores, additional arguments
(`_...`).

The package provides a number of macros to facilitate construction of
perturbation functions.  See [`@perturbn`](@ref `@perturbn`) and
[`@ivp`](@ref `@ivp`) in particular.  Thus for example the function
`p` above can be constructed so:
```
    p = @perturbn @lognormal(α,0.02) @lognormal(β,0.02)
```

## Initial value parameters

For certain types of parameters, one does not wish to apply the
perturbations at every lag. For example, parameters that specify the
initial conditions of the latent state process should have
perturbations applied only at lag 0. One can use the `lag` argument of
the `perturbations` function to accommodate these cases. For example,
suppose the perameters `α`, `β`, `γ`, and `δ` mentioned above are
regular parameters, but that `x₀` is a parameter that fixes the value
of the latent state at the zero-time. Then the following perturbations
function might be appropriate
```
    p(scale, lag; α, β, x₀, _...) = begin
        α=rand(LogNormal(log(α),scale*0.02))
        β=rand(LogNormal(log(β),scale*0.02))
        x₀ = (lag == 0) ? rand(LogNormal(log(x₀),scale*0.05)) : x₀
        (;α,β,x₀)
    end
```
Note that the perturbations are only applied to `x₀` at lag 0, i.e.,
at the zero-time.

Using the [`@perturbn`](@ref `@perturbn`) macro, the same function is
constructed via
```
    p = @perturbn(
            @lognormal(α,0.02),
            @lognormal(β,0.02),
            @ivp @lognormal(x₀,0.05)
        )
```

## Cooling schedule

The cooling schedule is specified by a function which, when furnished
a non-negative integer `n`, returns the fractional reduction of
perturbation scale (relative to that determined by the perturbation
kernel) that are to be applied in the `n`-th mif iteration.  The
function [`geometric_cooling`](@ref `geometric_cooling`) produces such
a function, but the user is free to specify alternative cooling
schedules.
"""
mif(
    object::ValidPompData;
    Np::Integer = 1,
    Nmif::Integer = 1,
    perturbations::Function,
    cooling::Function,
    trigger::Union{Real,Missing} = missing,
    target::Union{Real,Missing} = missing,
    params::P = coef(object),
    rinit::Union{Function,Nothing,Missing} = missing,
    rprocess::Union{PompPlugin,Nothing,Missing} = missing,
    logdmeasure::Union{Function,Nothing,Missing} = missing,
    kwargs...,
) where {P<:NamedTuple} = begin
    object = pomp(
        object;
        params=nothing,
        rinit, rprocess, logdmeasure,
        kwargs...,
    )
    trigger, target = proc_trig_targ(trigger, target)
    ll, trace = mif_internal(
        object, params,
        Nmif, Np,
        perturbations,
        cooling,
        trigger, target,
    )
    MifdPompObject(
        pfilter(object; params=trace[end], Np, trigger, target),
        Nmif,
        cooling,
        perturbations,
        trace, ll,
    )
end

mif(
    object::MifdPompObject;
    Np::Integer = object.pfobj.Np,
    Nmif::Integer = object.Nmif,
    perturbations::Function = object.perturbations,
    cooling::Function = object.cooling,
    trigger::Union{Real,Missing} = object.pfobj.trigger,
    target::Union{Real,Missing} = object.pfobj.target,
    kwargs...,
) = mif(
    pomp(object);
    Np, Nmif,
    perturbations, cooling,
    trigger, target,
    kwargs...,
)

mif(
    object::PfilterdPompObject;
    Np::Integer = object.Np,
    trigger::Union{Real,Missing} = object.trigger,
    target::Union{Real,Missing} = object.target,
    kwargs...,
) = mif(
    pomp(object);
    Np, trigger, target,
    kwargs...,
)

mif(_...) = error("Incorrect call to `mif`.")

mif_internal(
    object::PompObject{T,X,Y},
    params::P,
    Nmif::Integer,
    Np::Integer,
    perturbations::Function,
    cooling::Function,
    trigger::Union{Missing,Float64},
    target::Union{Missing,Float64},
) where {T,X,Y,P} = begin
    t0 = timezero(object)
    t = times(object)
    y = obs(object)
    params = fill(params,Np)
    ell = Array{LogLik}(undef,1,Np,1)
    cll = similar(ell,length(t))
    ess = Array{LogLik}(undef)
    pscale = cooling(0)
    perturbn!(params, perturbations, pscale, 0)
    x0 = rinit(object; params)
    xp = similar(x0, 1, Np, 1)
    xf = similar(x0, 1, Np, 1)
    resample = Array{Bool}(undef)
    perm = Array{Int}(undef,Np)
    work = similar(ell,Np)
    ll = similar(ell,Nmif)
    trace = Array{P}(undef,Nmif+1)
    trace[1] = param_mean(params)
    mif_calc!(
        trigger, target, Np,
        object, trace, ll, Nmif, x0, xp, xf, t0, t, y, params,
        ell, cll, ess, perturbations, cooling, pscale,
        perm, resample, work,
    )
    ll, trace
end

mif_calc!(
    trigger::Missing, target::Missing, Np::Integer, args...,
) = mif_loop!(args...)

mif_calc!(
    trigger::Float64, target::Float64, Np::Integer, args...,
) = begin
    w = ones(LogLik,Np)
    mif_loop!(args...,w,trigger,target,Np)
end

mif_loop!(
    object, trace, ll,
    Nmif, x0, xp, xf, t0, t, y, params,
    ell, cll, ess,
    perturbations, cooling, pscale,
    args...,
) = begin
    for i ∈ 1:Nmif
        t0 = timezero(object)
        if i > 1
            perturbn!(params, perturbations, pscale, 0)
            rinit!(object, x0; t0, params)
        end
        for j ∈ eachindex(t)
            mif_pfilt_step!(
                object,
                ell, @view(cll[j]), ess,
                xp, xf, params, x0,
                t0, @view(t[[j]]),
                @view(y[[j]]),
                args...,
            )
            if j < length(t) || i < Nmif
                perturbn!(params, perturbations, pscale, j)
            end
            t0 = t[j]
            x0 = view(xf,1,:,:)
        end
        ll[i] = sum(cll)
        trace[i+1] = param_mean(params)
        pscale = cooling(i)
    end
    nothing
end

mif_pfilt_step!(
    object::PompObject{T,X,Y},
    ell::AbstractArray{W,3},
    cll::AbstractArray{W,0},
    ess::AbstractArray{W,0},
    xp::AbstractArray{X,3},
    xf::AbstractArray{X,3},
    params::AbstractArray{P,1},
    x0::AbstractArray{X,2},
    t0::T,
    times::AbstractArray{T,1},
    y::AbstractArray{Y,1},
    perm::AbstractArray{I,1},
    args...,
) where {T,X,Y,W,P,I} = begin
    mif_advance_particles!(
        object, t0, times,
        x0, xp, y, ell, params,
    )
    pfilt_step_comps!(
        cll, ess,
        @view(ell[1,:,1]),
        perm,
        @view(xp[1,:,1]),
        @view(xf[1,:,1]),
        args...,
    )
    params .= params[perm]
    nothing
end

mif_advance_particles!(
    object, t0, times, x0, x, y, w, params,
) = begin
    rprocess!(object, x; x0, t0, times, params)
    logdmeasure!(object, w; times, y, x, params)
    nothing
end

## apply the perturbation kernel
perturbn!(
    params::Vector{P}, kernel::Function,
    scale::AbstractFloat, lag::Integer,
) where {P <: NamedTuple} = begin
    ptb = map(x -> kernel(scale, lag; x...), params)
    params .= P.(merge.(params,ptb))
    nothing
end

param_mean(
    params::Vector{P},
) where {P <: NamedTuple} = begin
    names = fieldnames(P)
    means = map(names) do n
        mean(getfield.(params,n))
    end
    (;zip(names,means)...)
end

pretty_string(object::MifdPompObject) = begin
    pretty_string(pomp(object)) *
        ", Nmif=$(object.Nmif)" *
        ", Np=$(object.pfobj.Np)" *
        ", logLik=$(round(logLik(object),digits=2))"
end

"""
    geometric_cooling(frac)

Returns a geometric cooling schedule under which the perturbations are
at a fraction `frac` of their original magnitude after 50 iterations.
"""
geometric_cooling(
    frac::AbstractFloat,
) = begin
    frac = Float64(frac)
    @assert 0 < frac ≤ 1 "`frac` must be ∈ (0,1]"
    speed = log(frac)/50
    n -> exp(speed*n)
end
