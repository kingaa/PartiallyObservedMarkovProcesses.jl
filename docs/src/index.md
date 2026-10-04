# PartiallyObservedMarkovProcesses.jl

The package is a Julia implementation of the [pomp package for R](https://kingaa.github.io/pomp/).

## Package Features

- [Implementation of POMP models](@ref)
- [Simulation](@ref)
- [Particle filter](@ref)
- [Iterated filtering](@ref)
- [Likelihood slices and profiles](@ref)
- [Plotting](@ref)
- [Trajectory matching](@ref)
- [Workhorses](@ref) (low-level interface to basic model components)
- [Helper functions](@ref)
- [Reproducibility tools](@ref)

## Function Documentation

### Implementation of POMP models

#### Basic constructor

```@docs
pomp
```

#### `rprocess` plugins

```@docs
euler
discrete_time
onestep
vectorfield
```

### Simulation

```@docs
simulate
simulate_array
```

### Particle filter

```@docs
pfilter
```

### Iterated filtering

```@docs
mif
@perturbn
geometric_cooling
hyperbolic_cooling
```

### Likelihood slices and profiles

The workflow follows R `pomp`: build a design with `slice_design` or
`profile_design`, evaluate it with `slice` (fixed-parameter particle
filters at every row) or `profile` (`mif` from every row, with the
profiled parameters held fixed, then fresh fixed-parameter particle
filters at the estimate), and summarize a profile with `mcap`.
`monitor` evaluates the unperturbed likelihood along a `mif` trace, as
a diagnostic of the fit.

```@docs
slice_design
profile_design
runif_design
sobol_design
pfilter_loglik
slice
profile
mcap
monitor
resampled
```

### Plotting

Plots are drawn with AlgebraOfGraphics, through a package extension:
load `AlgebraOfGraphics` and a Makie backend (e.g., `CairoMakie`) to use
them.  Each function returns the drawn figure; additional arguments are
passed to AlgebraOfGraphics' `draw`: `sliceplot` for the output of
`slice`, `mcapplot` for an `MCAP`, `traceplot` for `mif` results, and
`filterplot` for the effective sample size and conditional log
likelihood of `pfilter` or `mif` results (see the Reference page).

### Trajectory matching

```@docs
traj_match_objfun
```

### Workhorses

```@docs
rinit
rinit!
```

```@docs
rprocess
rprocess!
```

```@docs
rmeasure
```

```@docs
logdmeasure
logdmeasure!
```

```@docs
logdprior
logdprior!
```

```@docs
rprior
```

### Helper functions

```@docs
coef
obs
states
init_state
times
timezero
melt
logmeanexp
traces
```

### Reproducibility tools

```@docs
@freeze
@bake
```
