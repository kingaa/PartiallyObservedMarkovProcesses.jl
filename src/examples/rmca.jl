export rmca

import Distributions: Normal, LogNormal, logpdf

"""
    rmca(
        a = 1e-3, b = 1e-3, c = 1,
        r = 1, K = 1e4, m = 0.8,
        N0 = 8, P0 = 1,
        σ = 0.01,
        δt = 0.001, t₀ = 0.0,
        times=range(start=0.0,stop=500,step=0.1)
       )

`rmca` returns a *PompObject* containing simulated data.

## Parameters

- r: intrinsic growth rate of prey
- K: carrying capacity for prey
- a, c: predator foraging parameters
- b: predator yield
- N0, P0: initial conditions
- σ: measurement noise magnitude
- δt: Euler stepsize
- t₀: zero-time
- times: vector of observation times
"""
rmca = function(
    ;a = 1e-3, b = 1e-3, c = 1,
    r = 1, K = 1e4, m = 0.8,
    N0 = 8, P0 = 1,
    σ = 0.01,
    δt = 0.001, t₀ = 0.0,
    times=range(start=0,stop=500,step=0.1)
    )
    simulate(
        params=(
            ;a = a, b = b, c = c,
            r = r, K = K, m = m,
            N0 = N0, P0 = P0,
            σ = σ,
            δt = δt
        ),
        t0=t₀,
        times=times,
        rinit = function (;N0,P0,_...)
            (
                X=log(N0),
                Y=log(P0),
            )
        end,
        rprocess = function (;t,X,Y,a,b,c,r,K,m,δt,_...)
            N = exp(X)
            P = exp(Y)
            mu = [r*N*(1-N/K)*δt; c*P*N/(1+a*N)*δt; m*P*δt]
            dG = mu[1]*(1+randn())
            dS = mu[2]*(1+randn())
            dM = mu[3]*(1+randn())
            dG2 = mu[1]
            dS2 = mu[2]
            dM2 = mu[3]
            dX = (N > 0) ? (dG-dG2/N/2)/N-(dS-dS2/N/2)/N : -Inf
            dY = (P > 0) ? b*(dS-b*dS2/P/2)/P-(dM-dM2/P/2)/P : -Inf
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
