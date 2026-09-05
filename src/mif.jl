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
    "cooling function"
    cooling_schedule::Function
    "perturbation kernel"
    perturbation_kernel::Function
    "traces"
    trace::Array{P,1}
    "logLik trace"
    logLik::Array{W,1}
end

pomp(object::MifdPompObject) = pomp(object.pfobj)
logLik(object::MifdPompObject) = logLik(object.pfobj)
eff_sample_size(object::MifdPompObject) = eff_sample_size(object.pfobj)
cond_logLik(object::MifdPompObject) = cond_logLik(object.pfobj)

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
    mif(object; Np = 1, Nmif = 1, perturbation_kernel, cooling_schedule,
        trigger, target, params, rinit, rprocess, logdmeasure, kwargs...)

Iterated filtering.  In addition to the components needed for a
[`pfilter`](@ref `pfilter`) (i.e., `Np`, `trigger`, `target`), one
must specify a perturbation kernel, cooling schedule, and number of iterations.

## Arguments

- `object`: A `DataFrame`, `PompObject`, or vector of data.
- `Np`: number of particles to use in the filtering.
- `Nmif`: number of MIF iterations to perform.
- `perturbation_kernel`: a function that returns perturbed versions of
  some or all of the model parameters. See below for details.
- `cooling_schedule`: a function that specifies the MIF cooling
  schedule. See below for details.
- `trigger`, `target`: see [`pfilter`](@ref `pfilter`).
- `params`: `NamedTuple` of model parameters.
- `rinit`, `rprocess`, `logdmeasure`: necessary basic model components.
- `kwargs...`: other arguments are passed to [`pomp`](@ref `pomp`).

## Perturbation kernel

## Cooling schedule
"""
mif(
    object::ValidPompData;
    Np::Integer = 1,
    Nmif::Integer = 1,
    perturbation_kernel::Function,
    cooling_schedule::Function,
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
        perturbation_kernel,
        cooling_schedule,
        trigger, target,
    )
    MifdPompObject(
        pfilter(object; params=trace[end], Np, trigger, target),
        Nmif,
        cooling_schedule,
        perturbation_kernel,
        trace, ll,
    )
end

mif(
    object::MifdPompObject;
    Np::Integer = object.pfobj.Np,
    Nmif::Integer = object.Nmif,
    perturbation_kernel::Function = object.perturbation_kernel,
    cooling_schedule::Function = object.cooling_schedule,
    trigger::Union{Real,Missing} = object.pfobj.trigger,
    target::Union{Real,Missing} = object.pfobj.target,
    kwargs...,
) = mif(
    pomp(object);
    Np, Nmif,
    perturbation_kernel, cooling_schedule,
    trigger, target,
    kwargs...,
)

mif(_...) = error("Incorrect call to `mif`.")

mif_internal(
    object::PompObject{T,X,Y},
    params::P,
    Nmif::Integer,
    Np::Integer,
    perturbation_kernel::Function,
    cooling_schedule::Function,
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
    pscale = cooling_schedule(0)
    perturbn!(params, perturbation_kernel, pscale, 0)
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
        ell, cll, ess, perturbation_kernel, cooling_schedule, pscale,
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
    perturbation_kernel, cooling_schedule, pscale,
    args...,
) = begin
    for i ∈ 1:Nmif
        t0 = timezero(object)
        if i > 1
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
                perturbn!(params, perturbation_kernel, pscale, j)
            end
            t0 = t[j]
            x0 = view(xf,1,:,:)
        end
        ll[i] = sum(cll)
        trace[i+1] = param_mean(params)
        pscale = cooling_schedule(i)
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
    object::AbstractPompObject,
    t0::T,
    times::AbstractArray{T,1},
    x0::AbstractArray{X,2},
    x::AbstractArray{X,3},
    y::AbstractArray{Y,1},
    w::AbstractArray{W,3},
    params::AbstractArray{P,1},
) where {T,X,Y,W,P} = begin
    flexmap!(eachindex(params)) do j
        rprocess!(
            object, @view(x[:,[j],:]);
            t0, times,
            params=@view(params[[j]]),
            x0=@view(x0[[j],:]),
        )
        logdmeasure!(
            object, @view(w[:,[j],:]);
            times, y,
            params=@view(params[[j]]),
            x=@view(x[:,[j],:]),
        )
    end
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
