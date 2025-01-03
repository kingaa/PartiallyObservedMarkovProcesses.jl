using PartiallyObservedMarkovProcesses
using PartiallyObservedMarkovProcesses.Examples
using DataFrames
using Distributions
using Random
using RCall
using Test

@testset "speed trials" begin

    R"""
library(pomp)
library(tidyverse)
set.seed(599586410L)

P <- gompertz()
P |>
  as.data.frame() |>
  select(-X) -> dat

print(coef(P))

cat("pomp simulation times (Gompertz)\n")
P |>
  simulate(nsim=10000) |>
  system.time() |>
  getElement(3) |>
  replicate(n=3) |>
  print()

cat("pomp pfilter times (Gompertz)\n")
P |>
  pfilter(Np=10000,save.states="unweighted",filter.traj=TRUE) |>
  system.time() |>
  getElement(3) |>
  replicate(n=3) |>
  print()

cat("pomp likelihood estimate (Gompertz)\n")
P |>
  pfilter(Np=10000) |>
  logLik() |>
  round(digits=2) |>
  print()
""";

    @rget dat

    dat.time = Int64.(dat.time)

    theta = (K=1.0,r=0.1,σ=0.1,τ=0.1,X₀=1.0)

    P = pomp(
        dat,
        times=:time,
        t0=0,
        rinit=function(;X₀,_...)
            (X=X₀,)
        end,
        rprocess=discrete_time(
            function(;t,X,K,r,σ,_...)
                S = exp(-r)
                eps = if (σ > 0) exp(σ*randn()) else 1 end
                (X=K^(1-S)*X^S*eps,)
            end
        ),
        rmeasure=function(;X,τ,_...)
            d = LogNormal(log(X),τ)
            (Y=rand(d),)
        end,
        logdmeasure=function(;Y,X,τ,_...)
            d = LogNormal(log(X),τ)
            logpdf(d,Y)
        end
    )

    println("PartiallyObservedMarkovProcesses.jl simulation times (Gompertz)")
    simulate(P,params=theta,nsim=10000)
    @time simulate(P,params=theta,nsim=10000)
    @time simulate(P,params=theta,nsim=10000)
    @time simulate(P,params=theta,nsim=10000)

    println("PartiallyObservedMarkovProcesses.jl simulate scaling (Gompertz)")
    Q = simulate(P,params=theta,nsim=100);
    @time Q = simulate(P,params=theta,nsim=100)
    @time Q = simulate(P,params=theta,nsim=1000)
    @time Q = simulate(P,params=theta,nsim=10000)

    println("PartiallyObservedMarkovProcesses.jl simulate_array scaling (Gompertz)")
    simulate_array(P,params=theta,nsim=100);
    @time simulate_array(P,params=theta,nsim=100)
    @time simulate_array(P,params=theta,nsim=1000)
    @time simulate_array(P,params=theta,nsim=10000)

    println("PartiallyObservedMarkovProcesses.jl pfilter times (Gompertz)")
    pfilter(P,params=theta,Np=10000)
    @time pfilter(P,params=theta,Np=10000)
    @time pfilter(P,params=theta,Np=10000)
    @time pfilter(P,params=theta,Np=10000)

    println("PartiallyObservedMarkovProcesses.jl pfilter scaling (Gompertz)")
    @time Pf = pfilter(P,params=theta,Np=100)
    @time Pf = pfilter(P,params=theta,Np=1000)
    @time Pf = pfilter(P,params=theta,Np=10000)

    println(
        "PartiallyObservedMarkovProcesses.jl likelihood estimate (Gompertz): ",
        round(Pf.logLik,digits=2)
    )

end
