using PartiallyObservedMarkovProcesses
using PartiallyObservedMarkovProcesses.Examples
using DataFrames
using Distributions
using Random
using RCall
using BenchmarkTools

@info h1("speed trials")

@info h2("speed tests in R")

theta = (K=1.0,r=0.1,σ=0.1,τ=0.1,X0=1.0)

R"""
library(tidyverse,warn.conflicts=FALSE)
library(pomp,warn.conflicts=FALSE)
set.seed(599586410L)

$theta |>
  with(
    gompertz(K=K,r=r,sigma=σ,tau=τ,X_0=X0)
  ) -> P

P |>
  as.data.frame() |>
  select(-X) -> dat

cat("    pomp simulation times (Gompertz)\n")
P |>
  simulate(nsim=10000) |>
  system.time() |>
  getElement(3) |>
  replicate(n=3) |>
  print()

cat("    pomp pfilter times (Gompertz)\n")
P |>
  pfilter(Np=10000,save.states="filter",filter.traj=TRUE) |>
  system.time() |>
  getElement(3) |>
  replicate(n=3) |>
  print()

cat("    pomp likelihood estimate (Gompertz)\n")
P |>
  pfilter(Np=10000) |>
  logLik() -> ll

ll |>
  round(digits=2) |>
  print()
""";

@rget dat ll

dat.time = Int64.(dat.time)

P = pomp(
    dat,
    times=:time,
    t0=0,
    rinit=function(;X0,_...)
        (X=X0,)
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

@info h2("POMP.jl simulate scaling (Gompertz)")
@btime simulate($P,params=$theta,nsim=100)
@btime simulate($P,params=$theta,nsim=1000)
@btime simulate($P,params=$theta,nsim=10000)

@info h2("POMP.jl simulate_array scaling (Gompertz)")
@btime simulate_array($P,params=$theta,nsim=100)
@btime simulate_array($P,params=$theta,nsim=1000)
@btime simulate_array($P,params=$theta,nsim=10000)

@info h2("POMP.jl pfilter times (Gompertz)")
@btime pfilter($P,params=$theta,Np=10000)

@info h2("POMP.jl pfilter scaling (Gompertz)")
@btime pfilter($P,params=$theta,Np=100)
@btime pfilter($P,params=$theta,Np=1000)
@btime pfilter($P,params=$theta,Np=10000)
Pf = pfilter(P,params=theta,Np=10000)

@info h2("POMP.jl likelihood estimate (Gompertz): $(round(Pf.logLik,digits=2))")
@test abs(Pf.logLik-ll) < 1
