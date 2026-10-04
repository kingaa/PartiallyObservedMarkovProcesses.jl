import DataFrames: DataFrame, AbstractDataFrame, DataFrameRow, nrow, metadata, propertynames
import Random: default_rng

const LoglikSummary = NamedTuple{(:loglik,:se,:ess),NTuple{3,Float64}}

"""
    pfilter_loglik(object; Np, nreps = 1, params = coef(object), kwargs...)

Estimates the log likelihood at `params` from `nreps` replicate
[`pfilter`](@ref) runs, combined by [`logmeanexp`](@ref).  Returns
`(loglik, se, ess)`: the estimate, its jack-knife standard error, and
the effective sample size of the replicates.  Additional arguments are
passed to `pfilter`.
"""
pfilter_loglik(
    object::AbstractPompObject;
    Np::Integer,
    nreps::Integer = 1,
    params::NamedTuple = coef(object),
    kwargs...,
) = begin
    @assert nreps ≥ 1 "`nreps` must be positive"
    lls = [logLik(pfilter(object;Np,params,kwargs...)) for _ ∈ 1:nreps]
    summarize_loglik(lls)
end

pfilter_loglik(_...) = error("Incorrect call to `pfilter_loglik`.")

# A -Inf replicate represents a zero likelihood estimate and must be kept.
# Dropping it biases the likelihood average upward.
summarize_loglik(lls::AbstractVector{<:Real}) = begin
    @assert !any(isnan,lls) && !any(==(Inf),lls) "invalid NaN or +∞ log likelihood"
    if all(==(-Inf),lls)
        (loglik=-Inf,se=NaN,ess=0.0)
    elseif length(lls) == 1
        (loglik=Float64(lls[1]),se=NaN,ess=1.0)
    else
        est = logmeanexp(lls;se=true,ess=true)
        (loglik=est.est,se=est.se,ess=est.ess)
    end
end

# Override the parameters with the values in one design row.  The base is
# the model's stored parameters, updated by `params`; the row always wins.
# A model that stores parameters declares them: every design column except
# `slice` must then be one of them or be given in `params`.  A model that
# stores none declares nothing, and the design's columns are taken as given.
design_params(object::AbstractPompObject, row::DataFrameRow, params::NamedTuple = (;)) = begin
    stored = coef(object)
    base = merge(stored,params)
    for p ∈ keys(base)
        @assert p ∉ RESERVED_NAMES "parameter name `$p` is reserved: `slice`, `loglik`, `se`, and `ess` cannot be parameter names"
    end
    cols = filter(c -> c != :slice, propertynames(row))
    for c ∈ cols
        @assert c ∉ RESERVED_NAMES "design column `$c` is reserved for output"
    end
    if !isempty(stored)
        bad = filter(c -> c ∉ keys(base), cols)
        @assert isempty(bad) "design column(s) $(join(map(string,bad),", ")) are not parameters of the model"
    end
    merge(base,NamedTuple{Tuple(cols)}(Tuple(row[c] for c ∈ cols)))
end

with_loglik(design::AbstractDataFrame, res::AbstractVector{LoglikSummary}) = begin
    out = DataFrame(design;copycols=true)
    out.loglik = [r.loglik for r ∈ res]
    out.se = [r.se for r ∈ res]
    out.ess = [r.ess for r ∈ res]
    out
end

"""
    slice(object, design; Np, nreps = 1, params, kwargs...)

Likelihood slice.  Estimates the log likelihood by
[`pfilter_loglik`](@ref) at each row of `design` (typically from
[`slice_design`](@ref)), whose columns override the parameters of
`object`.  `params` supplies parameters the design does not contain,
on top of those stored in `object`; the design's values always win.  Returns `design` with `loglik`, `se`, and `ess` columns
appended, so these, and `slice`, cannot be parameter names.
Additional arguments are passed to `pfilter`.
"""
slice(
    object::AbstractPompObject,
    design::AbstractDataFrame;
    Np::Integer,
    nreps::Integer = 1,
    params::NamedTuple = (;),
    kwargs...,
) = begin
    n = nrow(design)
    @assert n > 0 "`design` has no rows"
    ps = [design_params(object,design[i,:],params) for i ∈ 1:n]
    res = Vector{LoglikSummary}(undef,n)
    flexmap!(1:n) do i
        res[i] = pfilter_loglik(object;Np,nreps,params=ps[i],kwargs...)
    end
    with_loglik(design,res)
end

slice(_...) = error("Incorrect call to `slice`.")

# The names of the parameters perturbed by a `mif` perturbations
# function, read from what it returns for each parameter set in `ps` at
# every lag `mif` uses (0 through `nlags`).  The state of the
# random-number generator is restored, so the fit is unchanged.
perturbed_names(
    perturbations::Function,
    ps::AbstractVector{<:NamedTuple},
    nlags::Integer,
) = begin
    rng = default_rng()
    saved = copy(rng)
    names = try
        union((keys(perturbations(1.0,lag;params...)) for params ∈ ps for lag ∈ 0:nlags)...)
    finally
        copy!(rng,saved)
    end
    Tuple(names)
end

# `perturbations`, checked at every call `mif` makes: any parameter it
# changes must be one of `perturbed`.  Probing cannot establish this
# for every function, e.g., one that depends on the cooling scale.
guarded(perturbations::Function, perturbed::Tuple) =
    function (scale, lag; params...)
        ptb = perturbations(scale,lag;params...)
        for p ∈ keys(ptb)
            @assert p ∈ perturbed || isequal(ptb[p],params[p]) "`perturbations` moved parameter `$p` during `mif` but did not return it when probed"
        end
        ptb
    end

"""
    profile(object, design; Nmif, Np, perturbations, cooling,
            nreps = 1, Np_eval = Np, profiled, params, kwargs...)

Profile likelihood.  From each row of `design` (typically from
[`profile_design`](@ref)), runs [`mif`](@ref), then estimates the log
likelihood at the resulting estimate by [`pfilter_loglik`](@ref) with
`Np_eval` particles and `nreps` replicates, using the model and the
resampling settings (`trigger`, `target`) that `mif` used.  Parameters not perturbed by
`perturbations` keep their values from the row.  Perturbing a profiled
parameter is an error: those recorded by `profile_design` and those
named in `profiled` are protected.  A design that carries no such
record (a plain data frame, or one widened with `hcat`, which drops
it) protects nothing unless `profiled` is given, and a warning says so;
`profiled = ()` declares that nothing is profiled.  `params` supplies
parameters the design does not contain, on top of those stored in
`object`; the design's values always win.  Additional arguments are
passed to `mif`.  `mif`'s default `avfun = geomean` fails for a
parameter that takes negative values, or that the model needs as an
integer (the geometric mean of a constant integer need not be exactly
that integer); pass `avfun = mean` then.  Returns a `DataFrame` of parameters
with `loglik`, `se`, and `ess` columns, suitable for [`mcap`](@ref); as
for [`slice`](@ref), these names are reserved.
"""
profile(
    object::AbstractPompObject,
    design::AbstractDataFrame;
    Nmif::Integer,
    Np::Integer,
    perturbations::Function,
    cooling::Function,
    nreps::Integer = 1,
    Np_eval::Integer = Np,
    profiled::Union{Nothing,Symbol,Tuple,AbstractVector} = nothing,
    params::NamedTuple = (;),
    kwargs...,
) = begin
    n = nrow(design)
    @assert n > 0 "`design` has no rows"
    ps = [design_params(object,design[i,:],params) for i ∈ 1:n]
    perturbed = perturbed_names(perturbations,ps,length(times(object)))
    # the profiled parameters: those the design records and those named
    md = metadata(design)
    record = haskey(md,"profiled") ? Symbol.(md["profiled"]) : nothing
    if isnothing(record) && isnothing(profiled)
        @warn "`profile`: the design does not record which parameters are profiled (a design from `profile_design` does, but `hcat` drops the record), so no parameter is protected against perturbation; pass `profiled = (...)`, or `profiled = ()` if none is."
    end
    named = profiled isa Symbol ? [profiled] : Symbol.(collect(something(profiled,())))
    protected = union(something(record,Symbol[]),named)
    for p ∈ protected
        @assert p ∈ keys(ps[1]) "profiled parameter `$p` is not a parameter"
        @assert p ∉ perturbed "profiled parameter `$p` must not be perturbed"
    end
    ests = Vector{NamedTuple}(undef,n)
    res = Vector{LoglikSummary}(undef,n)
    flexmap!(1:n) do i
        mf = mif(object;Nmif,Np,perturbations=guarded(perturbations,perturbed),cooling,params=ps[i],kwargs...)
        # unperturbed parameters keep their exact values: the swarm
        # average of a constant need not return that constant
        est = coef(mf)
        θ = merge(ps[i],NamedTuple{perturbed}(Tuple(getfield(est,p) for p ∈ perturbed)))
        ests[i] = θ
        # evaluate with the model and resampling settings (`trigger`,
        # `target`) that `mif` used, including any given in `kwargs`
        res[i] = pfilter_loglik(mf;Np=Np_eval,nreps,params=θ)
    end
    out = DataFrame(design;copycols=true)
    for p ∈ keys(ps[1])
        if p ∈ perturbed || p ∉ propertynames(out)
            out[!,p] = [getfield(ests[i],p) for i ∈ 1:n]
        end
    end
    with_loglik(out,res)
end

profile(_...) = error("Incorrect call to `profile`.")
