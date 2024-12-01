export rmca

import Distributions: Normal, LogNormal, logpdf

"""
    rmca(
        r = 1, K = 1e4, A = 1e3,
        b = 1e-3, c = 1, m = 0.8, V = 1000,
        N0 = 8, P0 = 1,
        σ = 0.01,
        δt = 0.01, t₀ = 0.0,
        times=range(start=0.0,stop=500,step=0.2)
       )

`rmca` returns a *PompObject* containing simulated data from a Rosenzweig-MacArthur model implemented as an Itô diffusion.
The process noise scales demographically (the system size is V) and the measurement error scales environmentally.

## Parameters

- r: intrinsic growth rate of prey
- K: carrying capacity for prey
- A, c: predator foraging parameters
- b: predator yield
- N₀, P₀: initial conditions
- V: system size
- σ: measurement noise magnitude
- δt: Euler stepsize
- t₀: zero-time
- times: vector of observation times

## Observables

- N: prey density
- P: predator density

## State variables

- X: ``X = log(N)``
- Y: ``Y = log(P)``
"""
rmca = function(
    ;r = 1, K = 1e4, A = 1e3,
    b = 1e-3, c = 1, m = 0.8,
    V = 100, σ = 0.01,
    N₀ = 3000, P₀ = 4,
    δt = 0.01, t₀ = 0.0,
    times=range(start=0,stop=500,step=0.2)
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
            mu = [r*N*δt; r*N*N/K*δt; c*P*N/(1+N/A)*δt; m*P*δt]
            eps = randn(4)/sqrt(V)
            dG = mu[1]+sqrt(mu[1])*eps[1]
            dC = mu[2]+sqrt(mu[2])*eps[2]
            dS = mu[3]+sqrt(mu[3])*eps[3]
            dM = mu[4]+sqrt(mu[4])*eps[4]
            dN = dG-dC-dS
            dP = b*dS-dM
            dN2 = (mu[1]+mu[2]+mu[3])/V/V
            dP2 = (b*b*mu[3]+mu[4])/V/V
            dX = (N > 0) ? (dN-dN2/N/2)/N : -Inf
            dY = (P > 0) ? (dP-dP2/P/2)/P : -Inf
            (
                t=t+δt,
                X=X+dX,
                Y=Y+dY
            )
        end,
        rmeasure = function (;X,Y,σ,_...)
            (
                N=rand(LogNormal(X,σ)),
                P=rand(LogNormal(Y,σ))
            )
        end,
        logdmeasure = function (;N,P,X,Y,σ,_...)
            logpdf(LogNormal(X,σ),N)+logpdf(LogNormal(Y,σ),P)
        end
    )[1]
end
