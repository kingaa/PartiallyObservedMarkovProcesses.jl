using PartiallyObservedMarkovProcesses
using PartiallyObservedMarkovProcesses.Examples
import PartiallyObservedMarkovProcesses as POMP
using Random
using RCall
using Test

println("- Gompertz model tests")

@testset verbose=true "Gompertz model" begin

    Random.seed!(1558102772)

    P = gompertz();
    @test isa(P,POMP.PompObject)
    println("    ",P)

    p1 = (r=4.5,K=210.0,σₚ=0.7,σₘ=0.1,X₀=150.0);

    x0 = rinit(P,nsim=5,params=p1);
    @test size(x0)==(1,5)
    @test keys(x0[1])==(:X,)

    x = rprocess(P,x0=x0,params=p1);
    @test size(x)==(27,1,5)
    @test keys(x[33])==(:X,)

    y = rmeasure(P,x=x,params=p1);
    @test size(y)==(27,1,5)
    @test keys(y[2])==(:pop,)

    y1 = rmeasure(P,x=x[3,:,:],times=1970,params=p1);
    @test size(y1)==(1,1,5)

    ell = logdmeasure(P,x=x,y=y,params=p1);
    @test size(ell)==(27,1,5,5)

    p = [p1; p1];
    x0 = rinit(P,params=p,nsim=3);
    x = rprocess(P,x0=x0,params=p);
    y = rmeasure(P,x=x[3:4,:,:],params=p,times=times(P)[3:4]);
    @test size(x0)==(2,3)
    @test size(x)==(27,2,3)
    @test size(y)==(2,2,3)

    Q = simulate(P,params=p1,nsim=3);
    Pf = pfilter(Q[1],Np=100);

    s = melt(Q,:parset,:rep);
    d = melt(P,rep=0);

    R"""
library(tidyverse,warn.conflicts=FALSE)
bind_rows($s,$d) |>
  mutate(data=if_else(rep==0,"data","simulation")) |>
  ggplot(aes(x=year,group=rep,color=factor(rep)))+
  geom_point(aes(y=pop,shape=data))+
  geom_line(data=$s,aes(y=X))+
  guides(color="none",size="none")+
  labs(shape="")+
  scale_y_sqrt()+
  theme_bw()
  """

    R"""ggsave(filename="gompertz-01.png",width=7,height=4)"""

end
