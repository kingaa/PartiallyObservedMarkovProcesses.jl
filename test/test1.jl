using POMP
using DataFrames
using CSV
using Distributions
using Random
using Test

dat = DataFrame(CSV.File("parus.csv",comment="#",delim=";"));

ri = function (;a,b,_...)
    d = NegativeBinomial(b,a)
    (x=rand(d), y=rand(d))
end

P = pomp(dat,times=:year,t0=0,params=(a=0.5,b=8));
@test rinit(P) broken=true

P = pomp(dat,times=:year,t0=0);
@test isa(coef(P),Nothing)

@test pomp(dat,times=:year,t0=1999) broken=true

P = pomp(dat,times=:year,t0=0,params=(a=0.5,b=8),rinit=ri);
@test isa(P,POMP.PompObject)

x = rinit(P)
@test isa(x,Vector)
@test length(x)==1
@test isa(x[1],NamedTuple)
@test length(x[1])==2
@test keys(x[1])==(:x,:y)

x = rinit(P,params=(b=9,a=0.1,))
@test isa(x,Vector)
@test length(x)==1
@test isa(x[1],NamedTuple)
@test length(x[1])==2
@test keys(x[1])==(:x,:y)

x = rinit(P,params=(b=9,a=0.2,),nsim=50)
@test isa(x,Vector)
@test length(x)==50
@test isa(x[10],NamedTuple)
@test length(x[10])==2
@test keys(x[10])==(:x,:y)

p = coef(P)
@test isa(p,NamedTuple)
@test length(p)==2
@test keys(p)==(:a,:b)

p = coef(P,:b,:a,:c)
@test isa(p,NamedTuple)
@test length(p)==2
@test keys(p)==(:b,:a)

pomp!(P,params=(a=0.1,b=200))
p = coef(P)
@test p==(a=0.1,b=200)

coef!(P,(a=13,c=22))
@test coef(P)==(a=13,c=22,b=200)
coef!(P,(a=13,c=22),reset=true)
@test coef(P)==(a=13,c=22)
coef!(P,reset=true)
@test isa(coef(P),Nothing)
