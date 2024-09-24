using POMP
using Random
using RCall
using Test

Random.seed!(1558102772)

P = gompertz()
@test isa(P,POMP.PompObject)
@test isnothing(states(P))
print(P)

x0 = rinit(P,nsim=5);
@test size(x0)==(1,5)
@test keys(x0[1])==(:x,)

x = rprocess(P,x0=x0);
@test size(x)==(27,1,5)
@test keys(x[33])==(:x,)

y = rmeasure(P,x=x);
@test size(y)==(27,1,5)
@test keys(y[2])==(:pop,)

@test_throws "with at least 3 dimensions" rmeasure(P,x=x[3,:,:],times=1970)

y1 = rmeasure(P,x=x[[3],:,:],times=1970);
@test size(y1)==(1,1,5)

p = [coef(P) for _ = 1:2];
x0 = rinit(P,params=p,nsim=3);
x = rprocess(P,x0=x0,params=p);
y = rmeasure(P,x=x[3:4,:,:],params=p,times=times(P)[3:4]);
@test size(x0)==(2,3)
@test size(x)==(27,2,3)
@test size(y)==(2,2,3)
@test_throws "x0-params size mismatch" rprocess(P,x0=x0,params=p[1])
@test_throws "x-times size mismatch" rmeasure(P,x=x[3:5,:,:],params=p,times=times(P)[3:4])
@test_throws "x-params size mismatch" rmeasure(P,x=x[3:4,:,:],params=p[2],times=times(P)[3:4])

coef!(P,(σₘ=0,σₚ=0))
coef(P)
melt(P);
X = rprocess(P,x0=rinit(P));
@test size(X)==(27,1,1)
d1 = melt(X,time=times(P));

@time Q = simulate(P;params=(r=4,σₘ=0,σₚ=0.7,x₀=0.1,K=210.0));
@test isa(Q,POMP.AbstractPompObject)
@test isa(Q,POMP.PompObject)
d2 = melt(Q);
@time Q = simulate(Q);
@time simulate!(Q);
@test isa(Q,POMP.PompObject)
d3 = melt(Q);
@test values(coef(Q,:σₘ,:σₚ)) == (0,0.7)
@time simulate!(Q;params=(r=4.5,σₘ=0.7,σₚ=0,x₀=0.1,K=210.0));
@time simulate!(Q);
@test values(coef(Q,:σₘ,:σₚ)) == (0.7,0)
@test isa(Q,POMP.PompObject)
d4 = melt(Q);

ell = dmeasure(Q);
@test isa(ell,Array{Float64,3})
@test size(ell)==(27,1,1)

x0 = rinit(Q,nsim=3);
@test size(x0)==(1,3)
x = rprocess(Q,x0=x0[:,1:2],times=timezero(Q).+[3,5,9]);
@test size(x)==(3,1,2)
y = rmeasure(Q,x=x[1:2,:,:],times=timezero(Q).+[3,5]);
@test size(y)==(2,1,2)
@test isa(obs(Q),Array{<:NamedTuple,3})
@test size(obs(Q))==(27,1,1)
@test_throws "replicate data not allowed" dmeasure(Q,y=y,x=x[1:2,:,:],times=timezero(Q).+[3,5]);
@test_throws "y-params size mismatch" dmeasure(Q,y=y,params=p,x=x[1:2,:,:],times=timezero(Q).+[3,5]);
@test_throws "y-times size mismatch" dmeasure(Q,y=y[[1],:,:],x=x[[1],:,:],times=timezero(Q).+[3,5]);
@test_throws "x-y size mismatch" dmeasure(Q,y=y[:,:,[1]],x=x[1:3,:,:],times=timezero(Q).+[3,5]);
@test_throws "at least 3 dimensions" dmeasure(Q,y=y[:,:,1],x=x[1:2,:,:],times=timezero(Q).+[3,5]);
@test_throws "at least 3 dimensions" dmeasure(Q,y=y,x=x[1:2,:,1],times=timezero(Q).+[3,5]);
ell = dmeasure(Q,y=y[:,:,[1]],x=x[1:2,:,:],times=timezero(Q).+[3,5]);
@test size(ell)==(2,1,2)

coef!(Q,(r=3,))
@test coef(Q,:r).r == 3
@test coef(Q).x₀ == 0.1

coef!(Q,(r=0.4,σₚ=0.2,x₀=0.1,K=1),reset=true)
@test_throws r"parameter .* undefined" dmeasure(Q)
coef!(Q,(σₘ=0,σₚ=0.2,x₀=0.1,K=1),reset=true)
@test_throws "parameter r undefined" rprocess(Q,x0=rinit(Q))
coef!(Q,(σₚ=0.2,x₀=0.1,K=1),reset=true)
@test_throws r"in `rmeasure`: parameter .* undefined" rmeasure(Q,x=states(Q))
coef!(Q,(σₚ=0.2,K=1),reset=true)
@test_throws r"in `rinit`: parameter .* undefined" rinit(Q)

@test_throws "basic component is undefined" simulate(P,rinit=nothing)
@test_throws "basic component is undefined" simulate(P,rprocess=nothing)
@test_throws "basic component is undefined" simulate(P,rmeasure=nothing)
@test_throws "in `simulate!`: in `rinit`: yikes!" simulate!(Q,rinit=function(;_...) error("yikes!") end)
@test_throws "in `simulate`: in `rprocess`: yikes!" simulate(P,rprocess=function(;_...) error("yikes!") end)
@test_throws "in `simulate`: in `rmeasure`: yikes!" simulate(P,rmeasure=function(;_...) error("yikes!") end)

R"""
library(tidyverse)
bind_rows(
    d1=$d1,
    d2=$d2,
    d3=$d3,
    d4=$d4,
    .id="sim"
  ) |>
  pivot_longer(-c(sim,time)) |>
  ggplot(aes(x=time,y=value,color=sim,linetype=name))+
  geom_point()+
  geom_line()+
  theme_bw()
"""

R"""ggsave(filename="gompertz-01.png",width=7,height=4)"""
