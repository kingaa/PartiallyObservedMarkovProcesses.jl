export rmca

import Distributions: Normal, LogNormal, logpdf

"""
    rmca(
         r = 1, K = 1e4, A = 1e3,
         b = 1e-3, c = 1, m = 0.8,
         V = 100, σ = 0.01,
         N₀ = 3000, P₀ = 4, t₀ = 0.0,
         δt = 0.01,
         times=range(start=0,stop=500,step=0.2)
       )


## Parameters

- r: intrinsic growth rate of prey
- K: carrying capacity for prey
- A: half-saturation prey density
- c: predator foraging rate
- b: predator yield (predators born per prey item killed)
- m: predator death rate
- V: system size
- σ: measurement noise magnitude
- N₀, P₀: initial densities
- t₀: zero-time
- δt: Euler stepsize
- times: vector of observation times

## Observables

- n: prey density
- p: predator density

## State variables

- ``X = \\log(N)``
- ``Y = \\log(P)``

## Details

`rmca` returns a *PompObject* containing simulated data from a Rosenzweig-MacArthur model implemented as an Itô diffusion.
Specifically, if ``N`` and ``P`` are prey and predator densities, respectively,
then ``dN = dG - dC - dS`` and ``dP = b dS - dM``, where
```math
\\begin{aligned}
dG &= r N dt + \\sqrt{\\frac{1}{V} r N} dW_1 \\\\
dC &= \\frac{r N^2}{K} dt + \\sqrt{\\frac{1}{V} \\frac{r N^2}{K}} dW_2 \\\\
dS &= \\frac{c N P}{1+N/A} dt + \\sqrt{\\frac{1}{V} \\frac{c N P}{1+N/A}} dW_3 \\\\
dM &= m P dt + \\sqrt{\\frac{1}{V} m P} dW_4 \\\\
\\end{aligned}
```
Here, the ``dW_i`` are increments of independent standard Wiener processes.
Thus, the process noise scales demographically.
Specifically, the system size, ``V``, converts the densities ``N``, ``P`` into numbers.
It controls the relative magnitude (coefficient of variation) of the demographic process noise.
Moreover, ``V`` determines a lower threshold on the population sizes, such that if ever ``N V < 1`` or ``P V < 1``, the population is taken to be extinct.
Otherwise, it plays no role in the dynamics.
The measurement error is assumed to scale environmentally:
```math
\\begin{aligned}
n &\\sim \\mathrm{LogNormal}(\\log{N},\\sigma) \\\\
p &\\sim \\mathrm{LogNormal}(\\log{P},\\sigma) \\\\
\\end{aligned}
```

Note that, in the limit ``V\\to\\infty``, the Itô diffusion becomes the ordinary differential equation
```math
\\begin{aligned}
\\frac{dN}{dt} &= r N \\left(1-\\frac{N}{K}\\right) - \\frac{c N P}{1+N/A} \\\\
\\frac{dP}{dt} &= \\frac{b c N P}{1+N/A} - m P \\\\
\\end{aligned}
```
which is the classical Rosenzweig-MacArthur model.

In this system, the predator is inviable unless ``R = \\frac{bcA}{m} > 1``.
Even if the predator is viable, the environment is too impoverished to support predators unless ``R>1+\\frac{A}{K}``.
If the environment is rich enough, and if moreover ``R>\\frac{1+\\frac{A}{K}}{1-\\frac{A}{K}}``, then the nontrivial equilibrium of the system is unstable.
For the default parameters, we have ``R = 1.25`` and ``\\frac{A}{K} = 0.1``, so the latter condition holds.
"""
rmca = function(
    ;r = 1, K = 1e4, A = 1e3,
    b = 1e-3, c = 1, m = 0.8,
    V = 100, σ = 0.01,
    N₀ = 3000, P₀ = 4, t₀ = 0.0,
    δt = 0.01, times=range(start=0,stop=500,step=0.2)
    )
    simulate(
        params=(
            ;r = r, K = K, A = A,
            b = b, c = c, m = m,
            N₀ = N₀, P₀ = P₀,
            V = V, σ = σ, δt = δt
        ),
        t0=t₀,
        times=times,
        rinit = function (;N₀,P₀,_...)
            (
                X=log(N₀),
                Y=log(P₀),
            )
        end,
        rprocess = function (;t,X,Y,b,c,r,K,A,m,V,δt,_...)
            N = exp(X)
            P = exp(Y)
            mu = [r*N; r*N*N/K; c*P*N/(1+N/A); m*P]
            dW = randn(4)*sqrt(δt/V)
            dG = mu*δt+sqrt.(mu).*dW
            dN = dG[1]-dG[2]-dG[3]
            dP = b*dG[3]-dG[4]
            dN2 = (mu[1]+mu[2]+mu[3])*δt/V
            dP2 = (b*b*mu[3]+mu[4])*δt/V
            dX = (N*V < 1) ? -Inf : (dN-dN2/N/2)/N
            dY = (P*V < 1) ? -Inf : (dP-dP2/P/2)/P
            (
                t=t+δt,
                X=X+dX,
                Y=Y+dY
            )
        end,
        rmeasure = function (;X,Y,σ,_...)
            (
                n=rand(LogNormal(X,σ)),
                p=rand(LogNormal(Y,σ))
            )
        end,
        logdmeasure = function (;n,p,X,Y,σ,_...)
            logpdf(LogNormal(X,σ),n)+logpdf(LogNormal(Y,σ),p)
        end
    )[1]
end
