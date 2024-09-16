using POMP
using Distributions
using Random
using DataFrames
using Test

dat = include("parus.jl");

rlni = function (;x₀,_...)
    (;x=x₀,)
end

rmeas = function (;x,σₘ,_...)
    d = LogNormal(log(x),σₘ)
    (;y=rand(d),)
end

P = pomp(
    dat,
    t0=0,    
    times=:year,
    params=(x₀=10,σₘ=0.1),
    rinit=rlni,
    rmeasure=rmeas
);
@test isa(P,POMP.PompObject)

x = rinit(P,nsim=5)
@test size(x)==(5,1,1)
@test keys(x[1])==(:x,)

y = rmeasure(P,state=x,time=3)
@test size(y)==(5,1,1)
@test keys(y[2])==(:y,)

x = rinit(P,params=[coef(P) for _ = 1:2],nsim=3)
y = rmeasure(P,state=x,params=[coef(P) for _ = 1:2],time=3)
@test size(x)==(3,2,1)
@test size(y)==(3,2,1)
@test keys(y[2])==(:y,)

