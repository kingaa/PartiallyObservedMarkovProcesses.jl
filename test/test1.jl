using POMP
using DataFrames
using Distributions
using Random
using RCall
using Test

Random.seed!(263260083)

rin = function(;x₀,_...)
    d = Poisson(x₀)
    (x=rand(d),)
end

rlin = function (;t,a,x,_...)
    d = Poisson(a*x)
    (t=t+1,x=rand(d))
end

rmeas = function (;x,k,_...)
    d = NegativeBinomial(k,k/(k+x))
    (y=rand(d),)
end

dmeas = function (;x,y,k,give_log,_...)
    d = NegativeBinomial(k,k/(k+x))
    if give_log
        logpdf(d,y)
    else
        pdf(d,y)
    end
end

@test_throws "must be no later than" pomp(times=[t for t in 0:20],t0=5)
@test_throws "times must be nondecreasing" pomp(t0=0,times=[-t for t in 0:7])
@test_throws "same elementary type" pomp(times=[t for t in 0:20],t0=0.0)

P = pomp(times=[t for t in 0:20],t0=0)
@test isa(P,POMP.PompObject)
@test timezero(P)==0
@test isa(times(P),Vector{<:Real})
@test length(times(P))==21
@test isnothing(obs(P))

p1 = (a=1.0,k=7.0,x₀=5.0);
p2 = (a=1.1,k=2.0,x₀=3.0);

x0 = rinit(P,params=p1);
@test isa(x0,Array{Tuple{},2})
@test size(x0)==(1,1)
@test x0[1]==()

pomp!(P,rinit=rin);
x0 = rinit(P,params=[p1;p2],nsim=7);
@test isa(x0,Array{<:NamedTuple,2})
@test size(x0)==(7,2)
@test keys(x0[8])==(:x,)
@test_throws r"parameter .* undefined" rinit(P,params=(a=1.0,k=7.0))

x = rprocess(P,x0=x0,params=[p1;p2]);
@test isa(x,Array{<:NamedTuple,3})
@test size(x)==(7,2,21)
@test x[:,:,21]==x0

pomp!(P,rprocess=rlin)
x = rprocess(P,x0=x0,params=[p1;p2]);
@test isa(x,Array{<:NamedTuple,3})
@test size(x)==(7,2,21)
@test x[:,:,1]==x0
@test x[:,:,21]!=x0
@test keys(x[17])==(:x,)
@test_throws r"parameter .* undefined" rprocess(P,x0=x0[1,:],params=(k=7.0,x₀=5.0))

y = rmeasure(P,x=x,params=[p1;p2]);
@test isa(y,Array{Tuple{},3})
@test size(y)==(7,2,21)
@test y[3]==()

P = pomp(P,rmeasure=rmeas)
y = rmeasure(P,x=x,params=[p1;p2]);
@test isa(y,Array{<:NamedTuple,3})
@test size(y)==(7,2,21)
@test isassigned(y,3)
@test keys(y[17])==(:y,)
@test_throws r"parameter .* undefined" rmeasure(P,x=x[1,:,:],params=(a=1.0,))

ell = dmeasure(P,x=x,y=y,params=[p1;p2]);
@test isa(ell,Array{Float64,4})
@test size(ell)==(size(y,1),size(x)...)
@test all(ell.==1.0)

ell = dmeasure(P,x=x,y=y,params=[p1;p2],give_log=true);
@test isa(ell,Array{Float64,4})
@test size(ell)==(size(y,1),size(x)...)
@test all(ell.==0.0)

pomp!(P,dmeasure=dmeas)
ell = dmeasure(P,x=x,y=y,params=[p1;p2],give_log=true);
@test isa(ell,Array{Float64,4})
@test size(ell)==(size(y,1),size(x)...)
@test all(ell.<=0)
@test_throws r"parameter .* undefined" dmeasure(P,y=y,x=x,params=(a=1.0,))

@test POMP.val_array("yes")==["yes"]
@test_throws "size mismatch" POMP.val_array(y,11,2)
@test size(POMP.val_array(3,1,1,1))==(1,1,1)

p = [p1;p2];
Q1 = simulate(P,params=p,nsim=3)
d1 = leftjoin(
    melt(Q1),
    melt(p,parset=1:2),
    on=:parset
);
Q2 = simulate(Q1);
d2 = leftjoin(
    melt(Q2),
    melt(p,parset=1:2),
    on=:parset
);
simulate!(Q2,nsim=4);
d3 = leftjoin(
    melt(p,parset=1:2),
    melt(Q2),
    on=:parset
);
simulate!(Q2,params=p2);
d4 = leftjoin(
    melt(Q2),
    melt(p2,parset=1),
    on=:parset
);
d = vcat(d1,d2,d3,d4,source=:sim);
R"""
library(tidyverse)
$d |>
  ggplot(aes(x=time,group=interaction(rep,sim),
    color=interaction(rep,sim)))+
  geom_line(aes(y=x))+
  geom_point(aes(y=y))+
  geom_line(aes(y=y),linetype=2)+
  guides(color="none")+
  facet_wrap(~a+k,ncol=1,scales="free_y",
    labeller=label_bquote(list(a==.(a),k==.(k))))+
  scale_y_continuous(transform="log1p")+
  theme_bw()
"""

R"""
ggsave(filename="test1-01.png",width=7,height=4)
"""

@test_throws "must be no later than" pomp(parus_data,times=:year,t0=1999)
@test_throws "times must be nondecreasing" pomp(sort(parus_data,:pop),times=:year,t0=1940)

pomp!(P,rinit=function(;_...) error("yikes!") end)
@test_throws "in `rinit`: yikes!" rinit(P,params=(x₀=3,))
pomp!(P,rprocess=function(;_...) error("yikes!") end)
@test_throws "in `rprocess!`: yikes!" rprocess(P,params=(x₀=3,),x0=x0)
pomp!(P,rmeasure=function(;_...) error("yikes!") end)
@test_throws "in `rmeasure`: yikes!" rmeasure(P,params=(a=1,),x=x)
pomp!(P,dmeasure=function(;_...) error("yikes!") end)
@test_throws "in `dmeasure!`: yikes!" dmeasure(P,params=(a=1,),x=x,y=y)

@test_throws "must be of the same length" pomp([(a=1,);(a=2,);(c=3,)];times=[1;2],t0=0)
