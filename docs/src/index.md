# PartiallyObservedMarkovProcesses.jl

The package is a Julia implementation of the [pomp package for R](https://kingaa.github.io/pomp/).

## Package Features

- [Implementation of POMP models](@ref)
- [Simulation](@ref)
- [Particle filter](@ref)
- [Trajectory matching](@ref)
- [Workhorses](@ref) (low-level interface to basic model components)
- [Helper functions](@ref)

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
```
