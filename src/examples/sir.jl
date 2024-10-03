export sir

"""
    sir()

`sir` is a *SimPompObject* containing simulated SIR data.

## Parameters
"""
sir = function(
    β = 0.5, γ = 0.25, N = 1000,
    ρ = 0.3, k = 5,
    S₀ = 0.6, I₀ = 0.01, R₀ = 0.4,
    δt = 0.1, t₀ = 0,
    times = [t for t ∈ range(0,90)],
    )
    simulate(
        params=(
            β=β,γ=γ,N=N,
            ρ=ρ,k=k,
            S₀=S₀,I₀=I₀,R₀=R₀,
            δt=δt,
        ),
        t0=t₀,
        times=times,
        rinit = function (;S₀,I₀,R₀,N,_...)
            m = N/(S₀+I₀+R₀)
            (
                S=round(Int64,m*S₀),
                I=round(Int64,m*I₀),
                R=round(Int64,m*R₀),
                C=0,
            )
        end,
        rprocess = function (;t,S,I,R,C,N,β,γ,δt,_...)
            infection = rand(Binomial(S,1-exp(-β*I/N*δt)))
            recovery = rand(Binomial(I,1-exp(-γ*δt)))
            (
                t=t+δt,
                S=S-infection,
                I=I+infection-recovery,
                R=R+infection,
                C=C+infection,
            )
        end,
        rmeasure = function (;ρ,k,C,_...)
            (
                reports=rand(NegativeBinomial(k,k/(k+ρ*C))),
            )
        end,
        logdmeasure = function (;reports,ρ,k,C,_...)
            logpdf(NegativeBinomial(k,k/(k+ρ*C)),reports)
        end
    )
end
