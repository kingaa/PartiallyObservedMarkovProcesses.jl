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

#### Basic constructor

```@docs
pomp
```

#### `rprocess` plugins

```@docs
euler
discrete_time
onestep
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
init_state
times
timezero
```

## Examples

```@setup
using POMP, RCall
R"""
options(tidyverse.quiet=TRUE)
library(tidyverse)
"""
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
  pivot_longer(-year) |>
  ggplot(aes(x=year,y=value))+
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
d = melt(Q,:parset,:rep)
R"""
svg("assets/figures/gompertz2.svg",width=7,height=5) #hide
$d |>
  pivot_longer(-c(year,rep,parset)) |>
  ggplot(aes(x=year,y=value,group=rep,color=factor(rep)))+
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

```@setup
using LatexSVG

texengine!(PDFLaTeX)
add_preamble!(
    "\\usepackage{tikz}",
    "\\definecolor{maize}{rgb}{1.0,0.796,0.020}",
    "\\definecolor{darkblue}{rgb}{0,0.153,0.298}",
    "\\definecolor{purple}{rgb}{0.518,0.082,0.765}",
    "\\definecolor{darkgreen}{rgb}{0,0.392,0}",
    "\\definecolor{royalblue}{rgb}{0.263,0.431,0.933}",
    "\\usetikzlibrary{arrows.meta,positioning,decorations,calc,math}",
)

sv = Lsvg"""
\resizebox{0.9\linewidth}{!}{
  \begin{tikzpicture}[scale=1]
    \tikzstyle{box}=[draw=black, text=black, fill=white, very thick, minimum size=3em]
    \tikzstyle{label}=[font=\Large]
    \tikzstyle{coordinate}=[inner sep=0pt,outer sep=0pt]
    \tikzstyle{flow}=[draw=black, very thick, >=stealth]
    \tikzstyle{modulate}=[draw=darkgreen, thick, >=Circle]
    \coordinate (origin) at (0,0);
    \node [box] (S) at ($(origin)+(1,-1)$) {${S}$};
    \node [box] (I) at ($(S)+(2,0)$) {${I}$};
    \node [box] (R) at ($(I)+(2,0)$) {${R}$};
    \coordinate (overR) at ($(R)+(0,1)$);
    \coordinate (midSI) at ($(S)!0.5!(I)$);
    \draw [flow,->] (S) -- (I);
    \draw [flow,->] (I) -- (R);
	%% \draw [flow,->] (R) -- (overR) -- (S |- overR) -- (S);
    \draw [modulate,->] (I.north west) .. controls ($(I)+(-0.5,0.5)$) and ($(midSI)+(0,1)$) .. (midSI);
  \end{tikzpicture}
}
"""
savesvg("assets/figures/sir_diagram.svg",sv,web_display=true)
```

![sir diagram](assets/figures/sir_diagram.svg)

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
