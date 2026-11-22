import DataFrames: DataFrame, metadata!, nrow
import Random: default_rng, seed!

"""
    hyperbolic_cooling(frac; start = 0)

Returns a hyperbolic cooling schedule under which the perturbations are
at a fraction `frac` of their original magnitude after 50 iterations:
the scale at iteration `n` is `a/(a+n)`, with `a = 50*frac/(1-frac)`.
Unlike R `pomp`'s, it does not cool within an iteration.  Continuing a
`mif` computation restarts the schedule; `start = k` resumes it after
`k` iterations.
"""
hyperbolic_cooling(
    frac::AbstractFloat;
    start::Integer = 0,
) = begin
    frac = Float64(frac)
    @assert 0 < frac ≤ 1 "`frac` must be ∈ (0,1]"
    @assert start ≥ 0 "`start` must be non-negative"
    if frac == 1
        n -> 1.0
    else
        a = 50*frac/(1-frac)
        n -> a/(a+n+start)
    end
end

"""
    resampled(object)

Returns the time series of indicators of whether resampling was performed in the particle filter computation stored in `object`.
"""
resampled(object::PfilterdPompObject) = object.resample
resampled(object::MifdPompObject) = resampled(object.pfobj)

"""
    monitor(runs; Np, seed, nreps = 1, every = 1)

Estimates the log likelihood by [`pfilter_loglik`](@ref) at the point
estimates recorded in the traces of one or more successive [`mif`](@ref)
computations, at every `every`-th iteration and the last.  Each run is
taken to continue the previous one (as `mif(previous; ...)` does), and
all must have the same parameters.  Unlike the
per-iteration log likelihoods in `traces(mf).logLik`, which are those of
the perturbed model, these are estimated by unperturbed particle
filters, each with the model and resampling settings of the
computation that recorded it.
They use their own random numbers, from `seed`, and leave those of the
session unchanged.  Returns a `DataFrame` of iterations, numbered as in
[`traces`](@ref), log likelihoods, and parameters; so parameters cannot
be named `iteration`, `loglik`, `se`, or `ess`.
"""
monitor(
    runs::AbstractVector{<:MifdPompObject};
    Np::Integer,
    seed::Integer,
    nreps::Integer = 1,
    every::Integer = 1,
) = begin
    @assert !isempty(runs) "no `mif` results given"
    @assert every ≥ 1 "`every` must be positive"
    pnames = keys(coef(runs[1]))
    for (k,m) ∈ enumerate(runs)
        d = symdiff(keys(coef(m)),pnames)
        @assert isempty(d) "the runs given to `monitor` must have the same parameters: run $k differs in $(join(map(p -> "`$p`",d),", "))"
    end
    for p ∈ pnames
        @assert p ∉ (:iteration,:loglik,:se,:ess) "parameter name `$p` is reserved: `iteration`, `loglik`, `se`, and `ess` cannot be parameter names"
    end
    # a continuation's first trace row is the perturbed start of that
    # run, not a completed iteration, so it is dropped
    points = NamedTuple[]
    models = AbstractPompObject[]
    iters = Int[]
    for (k,m) ∈ enumerate(runs)
        tr = traces(m)
        for j ∈ (k == 1 ? 1 : 2):nrow(tr)
            push!(points,NamedTuple{pnames}(Tuple(tr[j,p] for p ∈ pnames)))
            push!(models,m)
            push!(iters,isempty(iters) ? 1 : iters[end]+1)
        end
    end
    keep = [i for i ∈ eachindex(iters) if (iters[i]-1) % every == 0 || i == lastindex(iters)]
    rng = default_rng()
    saved = copy(rng)
    res = try
        seed!(seed)
        [pfilter_loglik(models[i];Np,nreps,params=points[i]) for i ∈ keep]
    finally
        copy!(rng,saved)
    end
    out = DataFrame(
        iteration=iters[keep],
        loglik=[r.loglik for r ∈ res],
        se=[r.se for r ∈ res],
        ess=[r.ess for r ∈ res],
    )
    for p ∈ pnames
        out[!,p] = [points[i][p] for i ∈ keep]
    end
    for (k,v) ∈ pairs((Np=Np,nreps=nreps,every=every,seed=seed))
        metadata!(out,string(k),v,style=:note)
    end
    out
end

monitor(run::MifdPompObject; kwargs...) = monitor([run]; kwargs...)

monitor(_...) = error("Incorrect call to `monitor`.")
