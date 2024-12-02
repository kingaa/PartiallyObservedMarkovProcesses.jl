# POMP.jl

The package is a Julia implementation of the [pomp package for R](https://kingaa.github.io/pomp/).

## Package Features

- [Implementation of POMP models](@ref)
- [Simulation](@ref)
- [Particle filter](@ref)
- [Workhorses](@ref) (low-level interface to basic model components)
- [Helper functions](@ref)
- [Examples](@ref)

## Function Documentation

```@setup
mkpath("assets/figures")
```

### Implementation of POMP models

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

## Examples

```@example
using POMP, RCall               #hide
R"""                            #hide
options(tidyverse.quiet=TRUE)   #hide
library(tidyverse)              #hide
"""                             #hide
nothing                         #hide
```

### The Gompertz model

```@docs
gompertz
```

View the Parus data:

```@example
using POMP, RCall
P = gompertz()
d = melt(P)
R"""
svg("assets/figures/gompertz1.svg",width=7,height=4) #hide
$d |>
  pivot_longer(-time) |>
  ggplot(aes(x=time,y=value))+
  geom_line()+
  facet_wrap(~name,scales="free_y",ncol=1)+
  labs(y="")+
  theme_bw() -> pl
print(pl)
dev.off() #hide
"""
nothing #hide
```

![gompertz_data](assets/figures/gompertz1.svg)

View a few representative simulations:

```@example
using POMP, RCall
P = gompertz()
Q = simulate(P;params=(r=4.5,K=210.0,σₚ=0.7,σₘ=0.1,X₀=150.0),nsim=5)
d = melt(Q,:rep,:parset)
R"""
svg("assets/figures/gompertz2.svg",width=7,height=5) #hide
$d |>
  pivot_longer(-c(time,rep,parset)) |>
  ggplot(aes(x=time,y=value,group=rep,color=factor(rep)))+
  geom_line()+
  geom_point()+
  facet_wrap(~name,scales="free_y",ncol=1)+
  labs(y="",color="replicate")+
  theme_bw() -> pl
print(pl)
dev.off() #hide
"""
nothing #hide
```

![gompertz_sims](assets/figures/gompertz2.svg)

### A simple SIR model

```@docs
sir
```

### The Rosenzweig-MacArthur model

```@docs
rmca
```

A sample simulation.

```@example
using POMP, RCall #hide
P = rmca(σ=0.1,times=range(0,400.0,step=1.0))
d = melt(P)
R"""
svg("assets/figures/rmca1.svg",width=7,height=6) #hide
$d |>
  pivot_longer(-time) |>
  ggplot(aes(x=time,y=value))+
  geom_path()+
  facet_wrap(~name,scales="free_y",ncol=1)+
  labs(y="")+
  theme_bw() -> pl
print(pl)
dev.off() #hide
"""
nothing #hide
```

![rmca_dynamics](assets/figures/rmca1.svg)
