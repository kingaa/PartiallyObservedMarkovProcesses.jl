import DataFrames: DataFrame, AbstractDataFrame, metadata!, metadata, nrow
import Sobol: SobolSeq, next!
import Random: AbstractRNG, default_rng

# column names that `slice` and `profile` use for their own output
const RESERVED_NAMES = (:slice,:loglik,:se,:ess)

# The values of a design column: integers stay integers (a model may need
# an integer parameter, e.g. a number of individuals); anything else
# becomes Float64.
design_values(v::Real) = design_values([v])
design_values(v) = begin
    v = collect(v)
    eltype(v) <: Integer ? v : collect(Float64,v)
end

# The values along a slice take the type of the parameter's value in
# `center` when that is an integer type and every value converts to it
# and back exactly (an integer parameter stays an integer); otherwise they
# become Float64.  The round trip, not a range check, decides: 2.0^63 is
# not ≤ typemax(Int) exactly, although the comparison in Float64 says so.
sliced_values(c::Real, vals) = begin
    v = vals isa Real ? [vals] : collect(vals)
    T = typeof(c)
    exact(x) = x isa Real && try T(x) == x catch; false end
    c isa Integer && all(exact,v) ? T.(v) : collect(Float64,v)
end

"""
    slice_design(center; slices...)

Design for a likelihood slice, as R `pomp`'s `slice_design`.  Each
keyword argument names a parameter of `center` and gives the values it
takes along its slice; the other parameters are held at `center`.
Returns a `DataFrame` with one column per parameter and a `slice`
column naming the parameter varied in each row.  A parameter whose
value in `center` is an integer stays an integer, held or sliced (when
the values along its slice are whole numbers), so a model that needs
an integer parameter gets one; other values become Float64.
"""
slice_design(center::NamedTuple; slices...) = begin
    @assert !isempty(slices) "at least one slice must be given"
    pnames = keys(center)
    @assert all(v -> v isa Real, values(center)) "`center` must be a `NamedTuple` of real numbers"
    for p ∈ pnames
        @assert p ∉ RESERVED_NAMES "parameter name `$p` is reserved: `slice`, `loglik`, `se`, and `ess` cannot be parameter names"
    end
    frames = map(collect(pairs(slices))) do (nm,vals)
        @assert nm ∈ pnames "variable `$nm` does not appear in `center`"
        v = sliced_values(getfield(center,nm),vals)
        n = length(v)
        df = DataFrame([p => design_values(fill(getfield(center,p),n)) for p ∈ pnames]...)
        df[!,nm] = v
        df[!,:slice] = fill(nm,n)
        df
    end
    reduce(vcat,frames)
end

slice_design(_...) = error("Incorrect call to `slice_design`.")

check_bounds(lower::NamedTuple, upper::NamedTuple) = begin
    ln = keys(lower)
    @assert Set(ln) == Set(keys(upper)) "names of `lower` and `upper` must match"
    @assert !isempty(ln) "`lower` and `upper` must name at least one variable"
    for k ∈ ln
        @assert lower[k] isa Real && upper[k] isa Real "`lower` and `upper` must be `NamedTuple`s of real numbers"
        @assert upper[k] ≥ lower[k] "upper values should be at least as large as lower ones (variable `$k`)"
    end
    ln
end

"""
    runif_design(lower, upper, nseq; rng = Random.default_rng())

`nseq` points drawn uniformly from the box with corners `lower` and
`upper`, as R `pomp`'s `runif_design`.  Returns a `DataFrame` with one
column per variable.
"""
runif_design(
    lower::NamedTuple,
    upper::NamedTuple,
    nseq::Integer;
    rng::AbstractRNG = default_rng(),
) = begin
    @assert nseq ≥ 0 "`nseq` must be non-negative"
    ln = check_bounds(lower,upper)
    DataFrame(
        [k => Float64(lower[k]) .+ (Float64(upper[k])-Float64(lower[k])).*rand(rng,nseq)
         for k ∈ ln]...
    )
end

runif_design(_...) = error("Incorrect call to `runif_design`.")

"""
    sobol_design(lower, upper, nseq)

`nseq` points of a Sobol' sequence spanning the box with corners
`lower` and `upper`, as R `pomp`'s `sobol_design`.  The sequence comes
from `Sobol.jl`, so the points differ from R's.
"""
sobol_design(
    lower::NamedTuple,
    upper::NamedTuple,
    nseq::Integer,
) = begin
    @assert nseq ≥ 0 "`nseq` must be non-negative"
    ln = check_bounds(lower,upper)
    d = length(ln)
    s = SobolSeq(d)
    x = zeros(Float64,d)
    m = Matrix{Float64}(undef,nseq,d)
    for i ∈ 1:nseq
        next!(s,x)
        m[i,:] .= x
    end
    DataFrame(
        [ln[j] => Float64(lower[ln[j]]) .+ (Float64(upper[ln[j]])-Float64(lower[ln[j]])).*m[:,j]
         for j ∈ 1:d]...
    )
end

sobol_design(_...) = error("Incorrect call to `sobol_design`.")

"""
    profile_design(; lower, upper, nprof, type = :runif, rng, profiled...)

Design for a profile likelihood, as R `pomp`'s `profile_design`.  The
keyword arguments `profiled...` give the values of the profiled
parameters; for each point of their grid, `nprof` starting values of
the other parameters are drawn from the box with corners `lower` and
`upper`, uniformly (`type = :runif`) or from a Sobol' sequence
(`type = :sobol`).  All values are Float64: a grid written with
integers, e.g. `K = 150:30:270`, is a grid of real values (`mif` averages
the parameters, and an integer column would not hold the average).  The profiled names are recorded in the metadata
under `"profiled"`, which [`profile`](@ref) checks.
"""
profile_design(;
    lower::NamedTuple,
    upper::NamedTuple,
    nprof::Integer,
    type::Symbol = :runif,
    rng::AbstractRNG = default_rng(),
    profiled...,
) = begin
    @assert !isempty(profiled) "at least one variable to profile over must be given"
    @assert nprof > 0 "`nprof` must be positive"
    @assert type ∈ (:runif,:sobol) "`type` must be `:runif` or `:sobol`"
    pv = collect(pairs(profiled))
    pnames = first.(pv)
    for k ∈ pnames
        @assert k ∉ keys(lower) "profiled variable `$k` must not appear in `lower`/`upper`"
    end
    for k ∈ (pnames...,keys(lower)...)
        @assert k ∉ RESERVED_NAMES "parameter name `$k` is reserved: `slice`, `loglik`, `se`, and `ess` cannot be parameter names"
    end
    grids = [v isa Real ? [Float64(v)] : collect(Float64,v) for (_,v) ∈ pv]
    rows = vec(collect(Iterators.product(map(eachindex,grids)...)))
    n = length(rows)
    x = DataFrame(
        [pnames[j] => [grids[j][r[j]] for r ∈ rows] for j ∈ eachindex(pnames)]...
    )
    y = if type == :runif
        runif_design(lower,upper,n*nprof;rng)
    else
        sobol_design(lower,upper,n*nprof)
    end
    out = hcat(repeat(x,inner=nprof),y)
    metadata!(out,"profiled",collect(pnames),style=:note)
    out
end

profile_design(_...) = error("Incorrect call to `profile_design`.")
