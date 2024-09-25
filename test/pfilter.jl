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

dmeas = function (;x,y,k,log,_...)
    d = NegativeBinomial(k,k/(k+x))
    if log
        logpdf(d,y)
    else
        pdf(d,y)
    end
end

p1 = (a=1.5,k=7.0,x₀=5.0);

P = simulate(
    t0=0,
    times=[i for i=1:20],
    params=p1,
    rinit=rin,
    rprocess=rlin,
    rmeasure=rmeas,
    dmeasure=dmeas
)

P = pomp(
    obs(P)[1,1,:],
    times=times(P),
    t0=0,
    rinit=rin,
    rprocess=rlin,
    rmeasure=rmeas,
    dmeasure=dmeas
)
@test isa(P,POMP.PompObject)

@time Q = pfilter(P,Np=1000,params=p1)
@time Q = pfilter(P,Np=1000,params=p1)
@time Q = pfilter(Q)
@test isa(Q,POMP.PfilterdPompObject)
@test size(states(Q))==(1000,1,20)
