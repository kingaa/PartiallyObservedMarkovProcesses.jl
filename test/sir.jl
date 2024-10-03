using POMP
using Random
using RCall
using Test

Random.seed!(1558102772)

P = sir()
@test isa(P,POMP.SimPompObject)
print(P)

@time Q = simulate(P,nsim=5)

d = melt(Q);

R"""
library(tidyverse)
$d |>
  select(-parset) |>
  pivot_longer(-c(rep,time)) |>
  ggplot(aes(x=time,y=value,group=rep,color=factor(rep)))+
  geom_line()+
  guides(color="none")+
  facet_wrap(~name,scales="free_y")+
  theme_bw()
ggsave(filename="sir-01.png",width=7,height=4)
"""
