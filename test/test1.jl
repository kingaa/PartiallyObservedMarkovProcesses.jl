using POMP
using DataFrames
using Distributions
using Random
using Test

dat = include("parus.jl");

rnb = function (;μ,k,_...)
    d = NegativeBinomial(k,k/(k+μ))
    (x=rand(d), y=rand(d))
end

P = pomp(dat,times=:year,t0=0,params=(μ=3,k=8));
@test_throws "basic component is undefined" rinit(P)

P = pomp(dat,times=:year,t0=0);
@test isa(coef(P),Nothing)

@test_throws "cannot be later than" pomp(dat,times=:year,t0=1999)
@test_throws "times must be nondecreasing" pomp(sort(dat,:pop),times=:year,t0=1940)

P = pomp(dat,times=:year,t0=0,params=(μ=10,k=0.1),rinit=rnb);
@test isa(P,POMP.PompObject)

x = rinit(P)
@test isa(x,Vector)
@test length(x)==1
@test isa(x[1],NamedTuple)
@test length(x[1])==2
@test keys(x[1])==(:x,:y)

x = rinit(P,params=(k=9,μ=0.1,))
@test isa(x,Vector)
@test length(x)==1
@test isa(x[1],NamedTuple)
@test length(x[1])==2
@test keys(x[1])==(:x,:y)

x = rinit(P,params=(μ=9,k=2,),nsim=50);
@test isa(x,Vector)
@test length(x)==50
@test isa(x[10],NamedTuple)
@test length(x[10])==2
@test keys(x[10])==(:x,:y)

p = coef(P)
@test isa(p,NamedTuple)
@test length(p)==2
@test keys(p)==(:μ,:k)

p = coef(P,:k,:μ,:r)
@test isa(p,NamedTuple)
@test length(p)==2
@test keys(p)==(:k,:μ)

pomp!(P,params=(k=0.1,μ=200))
p = coef(P)
@test p==(k=0.1,μ=200)

coef!(P,(μ=13,c=22))
@test coef(P)==(μ=13,c=22,k=0.1)
coef!(P,(μ=13,c=22),reset=true)
@test coef(P)==(μ=13,c=22)
coef!(P,reset=true)
@test isa(coef(P),Nothing)
