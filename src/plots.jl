# Plotting functions, implemented in the AlgebraOfGraphics extension
# (ext/PartiallyObservedMarkovProcessesAoGExt.jl): load
# `AlgebraOfGraphics` and a Makie backend (e.g., `CairoMakie`) to use them.

"""
    sliceplot(df; kwargs...)

Plots the output of [`slice`](@ref), one panel per sliced parameter,
with error bars from the `se` column.  Returns the figure drawn by
AlgebraOfGraphics; additional arguments are passed to `draw`.
Requires `AlgebraOfGraphics` and a Makie backend.
"""
function sliceplot end

"""
    mcapplot(m; kwargs...)

Plots an [`mcap`](@ref) result: the profile, the smoothed and quadratic
fits, the estimate, the confidence interval, and the cutoff.  Returns
the figure drawn by AlgebraOfGraphics; additional arguments are passed
to `draw`.  Requires `AlgebraOfGraphics` and a Makie backend.
"""
function mcapplot end

"""
    traceplot(mf; pars, monitor, kwargs...)

Plots the traces of one or a vector of [`mif`](@ref) computations, each
as a separate run: the log likelihood and the parameters `pars` (by
default, all those perturbed) against the iteration.  `monitor` adds a
panel from the output of [`monitor`](@ref), one data frame per run, each
computed for that run alone.  Returns the figure drawn by
AlgebraOfGraphics; additional arguments are passed to `draw`.  Requires
`AlgebraOfGraphics` and a Makie backend.
"""
function traceplot end

"""
    filterplot(x; kwargs...)

Plots the effective sample size and the conditional log likelihood
against time for the result of [`pfilter`](@ref) or [`mif`](@ref), or a
vector of them, one line per run.  For a `mif` result, these come from
the particle filter run at its final estimate.  Times at which the
conditional log likelihood is -∞ appear as gaps.  Returns the figure
drawn by AlgebraOfGraphics; additional arguments are passed to `draw`.
Requires `AlgebraOfGraphics` and a Makie backend.
"""
function filterplot end
