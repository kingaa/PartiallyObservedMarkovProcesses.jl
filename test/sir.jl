using POMP
using DataFrames
using Random
using RCall
using Test

@testset verbose=true "SIR model" begin

    Random.seed!(1558102772)

    R"""
library(tidyverse)
library(pomp)
P <- sir(
         gamma=0.25,mu=0,iota=0,rho=0.3,k=0.1,
         beta1=0.5,beta2=0.5,beta3=0.5,
         beta_sd=0,pop=10000,
         S_0=0.9,I_0=0.01,R_0=0.1,
         t0=0,times=seq(1,90,by=1)
         )
P |>
  as.data.frame() |>
  select(time,reports) -> dat

cat("pomp SIR parameters\n")
print(coef(P))

cat("pomp simulation times\n")
P |>
  simulate(nsim=100,format="a") |>
  system.time() |>
  getElement(3) |>
  replicate(n=3) |>
  print()

cat("pomp pfilter times\n")
P |>
  pfilter(Np=100,save.states="unweighted") |>
  system.time() |>
  getElement(3) |>
  replicate(n=3) |>
  print()

cat("pomp likelihood estimate\n")
P |>
  pfilter(Np=1000) |>
  logLik() -> ll
ll |>
  round(digits=2) |>
  print()
    """

    P = sir()
    @test isa(P,POMP.SimPompObject)
    print(P)

    println("POMP.jl simulation times")
    @time Q = simulate(P,nsim=100)
    @time Q = simulate(Q,nsim=100)
    simulate!(Q,accumvars=nothing)
    simulate!(Q)
    simulate!(Q)
    @time simulate!(Q,accumvars=(C=0,))
    @time simulate!(Q)
    @time simulate!(Q)
    @time simulate!(Q)
    @time simulate!(Q,nsim=5)

    d1 = melt(Q);

    R"""
library(tidyverse)
$d1 |>
  select(-parset) |>
  pivot_longer(-c(rep,time)) |>
  ggplot(aes(x=time,y=value,group=rep,color=factor(rep)))+
  geom_line()+
  guides(color="none")+
  facet_wrap(~name,scales="free_y")+
  theme_bw()
ggsave(filename="sir-01.png",width=7,height=4)
"""

    @rget dat

    println("POMP.jl SIR parameters")
    theta = coef(P)
    P = pomp(P);
    P.data = NamedTuple.(eachrow(select(dat,:reports)))

    println("POMP.jl pfilter times")
    @time Pf = pfilter(P,Np=100,params=theta)
    @time pfilter!(Pf)
    @time pfilter!(Pf)
    @time pfilter!(Pf)

    Pf = pfilter(Pf,Np=1000)
    println("POMP.jl likelihood estimate")
    println(round(Pf.logLik,digits=2))

    d2 = melt(Pf)

    R"""
library(tidyverse)
$d2 |>
  pivot_longer(-time) |>
  ggplot(aes(x=time,y=value))+
  geom_line()+
  guides(color="none")+
  facet_wrap(~name,scales="free_y",ncol=1)+
  theme_bw()
ggsave(filename="sir-02.png",width=7,height=4)
"""

end
