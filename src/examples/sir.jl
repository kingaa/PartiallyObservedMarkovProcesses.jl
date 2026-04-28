import Distributions: Binomial, NegativeBinomial, logpdf

"""
    sir(
        β = 0.5, γ = 0.25, N = 10000,
        ρ = 0.3, k = 10,
        S0 = 0.9, I0 = 0.01, R0 = 0.1,
        δt = 0.1, t0 = 0.0,
        times = range(start=1.0,stop=90,step=1.0)
       )

`sir` returns a *PompObject* containing simulated SIR data.

## Parameters

- β: transmission rate
- γ: recovery rate
- N: population size
- ρ: reporting rate
- k: overdispersion coefficient (negative binomial size parameter)
- S0, I0, R0: relative proportions of susceptible, infected, recovered (respectively) in the population at t=t0.
- δt: Euler stepsize
- t0: zero-time
- times: vector of observation times
"""
sir = function(
    ;β = 0.5, γ = 0.25, N = 10000,
    ρ = 0.3, k = 10,
    S0 = 0.9, I0 = 0.01, R0 = 0.1,
    δt = 0.1, t0 = 0.0,
    times = range(start=1.0,stop=90,step=1.0)
    )
    simulate(
        params = (
            β=β,γ=γ,N=N,
            ρ=ρ,k=k,
            S0=S0,I0=I0,R0=R0,
        ),
        t0 = t0,
        times = times,
        accumvars = (C=0,),
        rinit = function (;S0,I0,R0,N,_...)
            m = Float64(N)/(Float64(S0)+Float64(I0)+Float64(R0))
            (
                S=round(Int64,m*Float64(S0)),
                I=round(Int64,m*Float64(I0)),
                R=round(Int64,m*Float64(R0)),
                C=0,
            )
        end,
        rprocess = euler(
            function (;t,S,I,R,C,N,β,γ,dt,_...)
                infection = rand(Binomial(S,1-exp(-β*I/N*dt)))
                recovery = rand(Binomial(I,1-exp(-γ*dt)))
                (
                    S=S-infection,
                    I=I+infection-recovery,
                    R=R+recovery,
                    C=C+recovery,
                )
            end,
            dt = δt
        ),
        rmeasure = function (;ρ,k,C,_...)
            (
                reports=rand(NegativeBinomial(k,k/(k+ρ*C))),
            )
        end,
        logdmeasure = function (;reports,ρ,k,C,_...)
            logpdf(NegativeBinomial(k,k/(k+ρ*C)),reports)
        end
    )[1]
end
