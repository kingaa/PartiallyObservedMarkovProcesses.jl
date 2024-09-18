using POMP
using Distributions
using Random
using DataFrames
using RCall
using Test

Random.seed!(1558102772)

dat = include("parus.jl");

rlni = function (;x₀,_...)
    (;x=x₀,)
end

rmeas = function (;x,σₘ,_...)
    d = LogNormal(log(x),σₘ)
    (;y=rand(d),)
end

rproc = function (;t,x,σₚ,r,K,_...)
    s = exp(-r)
    d = LogNormal(s*log(x)+(1-s)*log(K),σₚ)
    (;t=t+1,x=rand(d),)
end

P = pomp(
    dat,
    t0=1960,
    times=:year,
    params=(x₀=0.1,σₘ=0.1,σₚ=0.2,r=0.1,K=1),
    rinit=rlni,
    rmeasure=rmeas,
    rprocess=rproc
);
@test isa(P,POMP.PompObject)

x = rinit(P,nsim=5)
@test size(x)==(5,1)
@test keys(x[1])==(:x,)

y = rmeasure(P,x=x,times=3)
@test size(y)==(5,1,1)
@test keys(y[2])==(:y,)

x = rinit(P,params=[coef(P) for _ = 1:2],nsim=3)
y = rmeasure(P,x=x,params=[coef(P) for _ = 1:2],times=3)
@test size(x)==(3,2)
@test size(y)==(3,2,1)
@test keys(y[2])==(:y,)

coef!(P,(r=0.2,σₘ=0,σₚ=0))
coef(P)
X = rprocess(P,x0=rinit(P));
d1 = hcat(DataFrame(t=times(P)),DataFrame(X))

Q = simulate(P;params=(r=0.3,σₘ=0,σₚ=0.2,x₀=0.1,K=1));
@test isa(Q,POMP.AbstractSimPompObject)
@test isa(Q,POMP.AbstractPompObject)
@test isa(Q,POMP.SimPompObject)
@test !isa(Q,POMP.PompObject)
d2 = hcat(DataFrame(t=times(Q)),DataFrame(states(Q)),DataFrame(obs(Q)))
Q = simulate(Q);
simulate!(Q)
@test isa(Q,POMP.SimPompObject)
d3 = hcat(DataFrame(t=times(Q)),DataFrame(states(Q)),DataFrame(obs(Q)))
@test values(coef(Q,:σₘ,:σₚ)) == (0,0.2)
simulate!(Q;params=(r=0.3,σₘ=0.1,σₚ=0,x₀=0.1,K=1));
@test values(coef(Q,:σₘ,:σₚ)) == (0.1,0)
@test isa(Q,POMP.SimPompObject)
d4 = hcat(DataFrame(t=times(Q)),DataFrame(states(Q)),DataFrame(obs(Q)))

x0 = rinit(Q,nsim=3);
@test size(x0)==(3,1)
x = rprocess(Q,x0=x0[1:2,:],times=timezero(Q).+[3,5,9]);
@test size(x)==(2,1,3)
y = rmeasure(Q,x=x[:,:,1:2],times=timezero(Q).+[3,5]);
@test size(y)==(2,1,2)
@test isa(obs(Q),Vector{<:NamedTuple})
@test length(obs(Q))==27

coef!(Q,(r=3,))
@test coef(Q,:r).r == 3
@test coef(Q).x₀ == 0.1

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
  pivot_longer(-c(sim,t)) |>
  mutate(t=unlist(t),value=unlist(value)) |>
  ggplot(aes(x=t,y=value,color=sim,linetype=name))+
  geom_point()+
  geom_line()+
  theme_bw()
"""
