import Distributions: Binomial, NegativeBinomial, logpdf

"""
    sir(
        β = 0.5, γ = 0.25, N = 10000,
        ρ = 0.3, k = 10,
        S₀ = 0.9, I₀ = 0.01, R₀ = 0.1,
        δt = 0.1, t₀ = 0.0,
        times = range(start=1.0,stop=90,step=1.0)
       )

`sir` returns a *PompObject* containing simulated SIR data.

## Parameters

- β: transmission rate
- γ: recovery rate
- N: population size
- ρ: reporting rate
- k: overdispersion coefficient (negative binomial size parameter)
- S₀, I₀, R₀: relative proportions of susceptible, infected, recovered (respectively) in the population at t=t₀.
- δt: Euler stepsize
- t₀: zero-time
- times: vector of observation times
"""
sir = function(
    ;β = 0.5, γ = 0.25, N = 10000,
    ρ = 0.3, k = 10,
    S₀ = 0.9, I₀ = 0.01, R₀ = 0.1,
    δt = 0.1, t₀ = 0.0,
    times = range(start=1.0,stop=90,step=1.0)
    )
    simulate(
        params = (
            β=β,γ=γ,N=N,
            ρ=ρ,k=k,
            S₀=S₀,I₀=I₀,R₀=R₀,
        ),
        t0 = t₀,
        times = times,
        accumvars = (C=0,),
        rinit = function (;S₀,I₀,R₀,N,_...)
            m = Float64(N)/(Float64(S₀)+Float64(I₀)+Float64(R₀))
            (
                S=round(Int64,m*Float64(S₀)),
                I=round(Int64,m*Float64(I₀)),
                R=round(Int64,m*Float64(R₀)),
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
