using POMP
using Random
using RCall
using Test

@testset "SIR model" begin

    Random.seed!(1558102772)

    P = sir()
    @test isa(P,POMP.SimPompObject)
    print(P)

    @time Q = simulate(P,nsim=5)
    @time Q = simulate(Q,nsim=5)
    @time simulate!(Q,accumvars=nothing)
    @time simulate!(Q)
    @time simulate!(Q)
    @time simulate!(Q,accumvars=(C=0,))
    @time simulate!(Q)
    @time simulate!(Q)

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

    dat = obs(P);
    theta = coef(P);
    P = pomp(P);
    P.data = vec(dat);

    @time Pf = pfilter(P,Np=1000,params=theta)
    @time pfilter!(Pf)
    @time pfilter!(Pf)
    @time pfilter!(Pf)

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
