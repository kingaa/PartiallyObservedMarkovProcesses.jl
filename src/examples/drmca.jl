import Distributions: LogNormal, logpdf
import DifferentialEquations: AutoTsit5, Rosenbrock23

"""
    drmca(
         r = 1, K = 1e4, A = 1e3,
         b = 1e-3, c = 1, m = 0.8,
         σ = 0.01,
         N₀ = 3000, P₀ = 4, t₀ = 0.0,
         times = range(start=0,stop=500,step=0.2),
         integrator = AutoTsit5(Rosenbrock23())
       )


## Parameters

- r: intrinsic growth rate of prey
- K: carrying capacity for prey
- A: half-saturation prey density
- c: predator foraging rate
- b: predator yield (predators born per prey item killed)
- m: predator death rate
- σ: measurement noise magnitude
- N₀, P₀: initial densities
- t₀: zero-time
- times: vector of observation times
- integrator: integration algorithm.

## Observables

- n: prey density
- p: predator density

## State variables

- ``X = \\log(N)``
- ``Y = \\log(P)``

## Details

`drmca` returns a *PompObject* containing simulated data from a deterministic Rosenzweig-MacArthur model implemented as a vectorfield.
Specifically, the model is the classical Rosenzweig-MacArthur model
```math
\\begin{aligned}
\\frac{dN}{dt} &= r N \\left(1-\\frac{N}{K}\\right) - \\frac{c N P}{1+N/A} \\\\
\\frac{dP}{dt} &= \\frac{b c N P}{1+N/A} - m P \\\\
\\end{aligned}
```

In this system, the predator is inviable unless ``R = \\frac{bcA}{m} > 1``.
Even if the predator is viable, the environment is too impoverished to support predators unless ``R>1+\\frac{A}{K}``.
If the environment is rich enough, and if moreover ``R>\\frac{1+\\frac{A}{K}}{1-\\frac{A}{K}}``, then the nontrivial equilibrium of the system is unstable.
For the default parameters, we have ``R = 1.25`` and ``\\frac{A}{K} = 0.1``, so the latter condition holds.
"""
drmca = function(
    ;r = 1, K = 1e4, A = 1e3,
    b = 1e-3, c = 1, m = 0.8,
    σ = 0.01,
    N₀ = 3000, P₀ = 4, t₀ = 0.0,
    times = range(start=0,stop=500,step=0.2),
    integrator = AutoTsit5(Rosenbrock23())
    )
    simulate(
        params=(
            ;r = r, K = K, A = A,
            b = b, c = c, m = m,
            N₀ = N₀, P₀ = P₀,
            σ = σ
        ),
        t0=t₀,
        times=times,
        rinit = function (;N₀,P₀,_...)
            (
                X=log(N₀),
                Y=log(P₀),
            )
        end,
        rprocess = vectorfield(
            function (X,Y;t,b,c,r,K,A,m,_...)
                N = exp(X)
                P = exp(Y)
                mu = [r*N, r*N*N/K, c*P*N/(1+N/A), m*P]
                dN = mu[1]-mu[2]-mu[3]
                dP = b*mu[3]-mu[4]
                [dN/N, dP/P]
            end,
            integrator
        ),
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
