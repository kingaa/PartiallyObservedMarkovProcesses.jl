import Statistics: mean
import DataFrames: insertcols!

"""
    MifdPompObject

Created by a call to [`mif`](@ref `mif`), this struct holds the
results of the particle-filter computation.  In particular, one can
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

pfilter(
    object::MifdPompObject;
    kwargs...,
) = pfilter(object.pfobj; kwargs...)

traces(object::MifdPompObject) = begin
    insertcols!(melt(object.trace, :iteration), :logLik => [object.logLik...,missing])
end

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
    t0 = timezero(object)
    t = times(object)
    y = val_array(obs(object),length(t),1)
    trace = Array{P}(undef,Nmif+1)
    params = fill(params,Np)
    trace[1] = param_mean(params)
    ell = Array{LogLik}(undef,length(t),Np,1,1)
    cll = similar(ell,length(t))
    ess = similar(cll)
    ll = similar(ell,Nmif)
    work = similar(ell,Np)
    perm = Array{Int}(undef,Np)
    pscale = cooling_schedule(0)
    perturbn!(params, perturbation_kernel, pscale, 0)
    x0 = POMP.rinit(object; params)
    xp = similar(x0, 1, Np, 1)
    xf = similar(x0, 1, Np, 1)
    resample = Array{Bool}(undef)
    for i ∈ 1:Nmif
        t0 = timezero(object)
        if i > 1
            rinit!(object, x0; t0, params)
        end
        for j ∈ eachindex(t)
            rprocess!(object, xp; params, x0, t0=t0, times=@view(t[[j]]))
            logdmeasure!(object, @view(ell[[j],:,:,:]); times=@view(t[[j]]), y = @view(y[[j],:,:]), x=xp, params)
            pfilt_step_comps!(
                @view(cll[j]), @view(ess[j]),
                @view(ell[j,:,1,1]),
                perm,
                @view(xp[1,:,1]),
                @view(xf[1,:,1]),
                resample,
                work,
            )
            params .= params[perm]
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
    Np, Nmif, perturbation_kernel, cooling_schedule,
    trigger, target, kwargs...,
)

mif(_...) = error("Incorrect call to `mif`.")

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
