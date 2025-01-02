using PartiallyObservedMarkovProcesses
using PartiallyObservedMarkovProcesses.Examples
using DataFrames
using Distributions
using Random
using RCall
using Test
# using Revise
using Debugger

import PartiallyObservedMarkovProcesses as POMP 

@testset "speed trials 2" begin

    NSIM = 100_000
    # NSIM = 1_000_000
    # NSIM = 1_000 
    NP = 10_000
    theta = (K=1.0,r=0.1,σ=0.1,τ=0.1,X₀=1.0)

    R"""
library(pomp)
library(tidyverse)
set.seed(599586410L)

P <- gompertz()
P |>
  as.data.frame() |>
  select(-X) -> dat

print(coef(P))
""";

    @rget dat

    dat.time = Int64.(dat.time)

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
            # d = LogNormal(log(X),τ)
            # (Y=rand(d),)
            (Y=1 ,)
        end,
        logdmeasure=function(;Y,X,τ,_...)
            # d = LogNormal(log(X),τ)
            # logpdf(d,Y)
            1 
        end
    )

    # POMP.setusethreads!(false)
    # @run simulate(P,params=theta,nsim=NSIM)
     
    POMP.setusethreads!(false)
    println("PartiallyObservedMarkovProcesses.jl simulation times (Gompertz)")
    @time simulate(P,params=theta,nsim=NSIM)

    println("PartiallyObservedMarkovProcesses.jl pfilter times (Gompertz)")
    pfilter(P,params=theta,Np=NSIM)
    @time pfilter(P,params=theta,Np=NP)


    POMP.setusethreads!(true)
    println("Multi-threaded PartiallyObservedMarkovProcesses.jl simulation times (Gompertz)")
    simulate(P,params=theta,nsim=NSIM)
    @time simulate(P,params=theta,nsim=NSIM)

    println("Multi-threaded PartiallyObservedMarkovProcesses.jl pfilter times (Gompertz)")
    pfilter(P,params=theta,Np=NSIM)
    @time pfilter(P,params=theta,Np=NP)

end
