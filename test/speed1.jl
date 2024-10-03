using POMP
using DataFrames
using Distributions
using Random
using RCall
using Test

R"""
library(pomp)
library(tidyverse)
set.seed(599586410L)

P <- gompertz()
P |>
  as.data.frame() |>
  select(-X) -> dat

print(coef(P))

cat("pomp simulation times\n")
P |>
  simulate(nsim=10000) |>
  system.time() |>
  getElement(3) |>
  replicate(n=3) |>
  print()

cat("pomp pfilter times\n")
P |>
  pfilter(Np=10000) |>
  system.time() |>
  getElement(3) |>
  replicate(n=3) |>
  print()

cat("pomp likelihood estimate\n")
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
    rprocess=function(;t,X,K,r,σ,_...)
        S = exp(-r)
        eps = if (σ > 0) exp(σ*randn()) else 1 end
        (t=t+1,X=K^(1-S)*X^S*eps,)
    end,
    rmeasure=function(;X,τ,_...)
        d = LogNormal(log(X),τ)
        (Y=rand(d),)
    end,
    logdmeasure=function(;Y,X,τ,_...)
        d = LogNormal(log(X),τ)
        logpdf(d,Y)
    end
)

println("POMP.jl simulation times")
@time simulate(P,params=theta,nsim=10000)
@time simulate(P,params=theta,nsim=10000)
@time Q = simulate(P,params=theta,nsim=10000)
@time simulate!(Q)
@time simulate!(Q)
@time simulate!(Q)

println("POMP.jl pfilter times")
@time Pf = pfilter(P,params=theta,Np=10000)
@time Pf = pfilter(P,params=theta,Np=10000)
@time Pf = pfilter(P,params=theta,Np=10000)
@time pfilter!(Pf)
@time pfilter!(Pf)
@time pfilter!(Pf)

println("POMP.jl likelihood estimate")
println(round(Pf.logLik,digits=2))
