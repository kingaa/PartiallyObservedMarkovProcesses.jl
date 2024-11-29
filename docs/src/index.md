# POMP.jl

Summary.

## Package Features

- Implement POMP models.
- Simulate models with process noise and measurement error.
- Particle filters.

## Function Documentation

### Basic constructor

```@docs
pomp
```

### Simulation

```@docs
simulate
```

### Particle filter

```@docs
pfilter
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

### Helper functions

```@docs
coef
obs
states
times
timezero
```

### Examples

```@docs
gompertz
sir
rmca
```
