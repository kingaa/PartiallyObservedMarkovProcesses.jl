using PartiallyObservedMarkovProcesses
using PartiallyObservedMarkovProcesses.Examples
import PartiallyObservedMarkovProcesses as POMP
using Random
using RCall
using Test
using BenchmarkTools

@info h1("SIR model tests")

@testset verbose=true "SIR model" begin

    Random.seed!(1558102772)

    theta = (γ=0.25,ρ=0.3,k=10,β=0.5,N=10000,S0=0.9,I0=0.01,R0=0.1);
    @info "SIR parameters: $theta"

    @info h2("tests with R pomp")

    R"""
library(tidyverse,warn.conflicts=FALSE)
library(pomp,warn.conflicts=FALSE)

set.seed(1558102772)

theta <- unlist($theta)

simulate(
  t0=0,
  times=seq(1,90,by=1),
  rinit=Csnippet(r"{
    double m = N/(S0+I0+R0);
    S = nearbyint(m*S0);
    I = nearbyint(m*I0);
    R = nearbyint(m*R0);
    C = 0;}"
  ),
  rprocess=euler(
    Csnippet(r"{
    double inf = rbinom(S,1-exp(-β*I/N*dt));
    double rcv = rbinom(I,1-exp(-γ*dt));
    S -= inf;
    I += inf - rcv;
    R += rcv;
    C += rcv;}"
    ),
    delta.t=0.1
  ),
  rmeasure=Csnippet(
    r"{reports = rnbinom(k,k/(k+ρ*C));}"
  ),
  dmeasure=Csnippet(
    r"{
      lik = dnbinom(reports,k,k/(k+ρ*C),give_log);
    }"
  ),
  obsnames=c("reports"),
  statenames=c("S","I","R","C"),
  paramnames=c(
    "γ","β","k","ρ","N",
    "S0","I0","R0"
  ),
  accumvars="C",
  params = theta
) -> P

P |>
  as.data.frame() |>
  select(time,reports) -> dat

cat("    pomp simulation times (SIR):\n")
P |>
  simulate(nsim=1000,format="p") |>
  system.time() |>
  getElement(3) |>
  replicate(n=3) |>
  print()

pp <- parmat(theta,2)
pp[c("γ","β"),2] <- c(0.5,0.25)
P |>
  simulate(nsim=1000,format="a",params=pp) |>
  system.time() |>
  getElement(3) |>
  replicate(n=3) |>
  print()

cat("    pomp pfilter times (SIR):\n")
P |>
  pfilter(Np=1000,save.states="filter",filter.traj=TRUE) |>
  system.time() |>
  getElement(3) |>
  replicate(n=3) |>
  print()

P |>
  pfilter(Np=1000) |>
  logLik() -> ll
    """;

    @rget ll
    @info h2("pomp likelihood estimate (SIR): $(round(ll,digits=2))")

    P = sir();
    @test P isa POMP.PompObject
    @test occursin(r"PompObject with 90 observations",sprint(show,P))

    @info h2("POMP.jl simulation times (SIR)")
    Q = simulate(P,nsim=1000,params=theta);
    @test typeof(Q[1])!=typeof(P)
    @btime simulate($(Q[1]),nsim=1000);

    Q = simulate(P,nsim=5);
    @test typeof(Q[1])==typeof(P)

    d = simulate_array(
        P,
        nsim=5,
        params=[
            (γ=0.25,ρ=0.3,k=10,β=0.5,N=10000,S0=0.9,I0=0.01,R0=0.1);
            (γ=0.5,ρ=0.3,k=10,β=0.25,N=10000,S0=0.9,I0=0.01,R0=0.1);
        ]
    );
    @test size(d)==(90,2,5)
    @test keys(d[1])==(:time,:reports,:S,:I,:R,:C)
    @test_throws "invalid Array dimensions" simulate_array(P,nsim=-1)

    @info h2("POMP.jl simulate_array times (SIR)")
    @btime simulate_array(
        $P,nsim=5,
        params=[
            (γ=0.25,ρ=0.3,k=10,β=0.5,N=10000,S0=0.9,I0=0.01,R0=0.1);
            (γ=0.5,ρ=0.3,k=10,β=0.25,N=10000,S0=0.9,I0=0.01,R0=0.1);
        ]
    );

    d1 = melt(Q,:parset,:rep);

    R"""
library(tidyverse,warn.conflicts=FALSE)
$d1 |>
  select(-parset) |>
  pivot_longer(-c(rep,time)) |>
  ggplot(aes(x=time,y=value,group=rep,color=factor(rep)))+
  geom_line()+
  guides(color="none")+
  facet_wrap(~name,scales="free_y")+
  theme_bw()
"""

    R"""ggsave(filename="sir-01.png",width=7,height=4)"""

    @rget dat
    P = pomp(
        dat,
        times=:time,
        t0=0.0,
        rinit=P.rinit,
        rprocess=P.rprocess,
        rmeasure=P.rmeasure,
        logdmeasure=P.logdmeasure,
        params=coef(P),
        accumvars=(C=0,)
    )

    @info h2("POMP.jl pfilter times (SIR)")
    Pf = pfilter(P,Np=1000,params=theta);
    @btime pfilter($P,Np=1000,params=$theta);

    Pf = pfilter(Pf,Np=1000);
    @info h2("POMP.jl likelihood estimate (SIR): $(round(Pf.logLik,digits=2))")
    @test abs(Pf.logLik-ll) < 1.0

    d2 = melt(Pf)

    R"""
library(tidyverse,warn.conflicts=FALSE)
$d2 |>
  pivot_longer(-time) |>
  ggplot(aes(x=time,y=value))+
  geom_line()+
  guides(color="none")+
  facet_wrap(~name,scales="free_y",ncol=1)+
  theme_bw()
"""

    R"""ggsave(filename="sir-02.png",width=7,height=4)"""

end
