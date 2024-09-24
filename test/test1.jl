using POMP
using DataFrames
using Distributions
using Random
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
    (pop=rand(d),)
end

dmeas = function (;x,pop,k,_...)
    d = NegativeBinomial(k,k/(k+x))
    logpdf(d,pop)
end

P = pomp(parus_data,times=:year,t0=1960)
@test isa(P,POMP.PompObject)

P = pomp(parus_data,times=:year,t0=1960)
@test isa(P,POMP.PompObject)
@test timezero(P)==1960
@test isa(times(P),Vector{<:Real})
@test length(times(P))==27
y=obs(P);
@test isa(y,Array{<:NamedTuple,3})
@test size(y)==(27,1,1)
@test length(y[1])==1
@test keys(y[1])==(:pop,)
@test_throws "not defined" obs!(P,y)

@test_throws "cannot be later than" pomp(parus_data,times=:year,t0=1999)
@test_throws "times must be nondecreasing" pomp(sort(parus_data,:pop),times=:year,t0=1940)

@test_throws "basic component is undefined" rinit(P,params=(a=1,k=3,x0=5))
P = pomp(parus_data,times=:year,t0=1960,rinit=rin);
@test isa(P,POMP.PompObject)

x0 = rinit(P,params=[(x₀=9,),(x₀=3,)],nsim=7);
@test isa(x0,Array)
@test size(x0)==(2,7)
@test isa(x0[10],NamedTuple)
@test length(x0[11])==1
@test keys(x0[11])==(:x,)

x0 = rinit(P,params=(a=1,k=3,x₀=5),nsim=5);
@test isa(x0,Array)
@test length(x0)==5
@test size(x0)==(1,5)
@test isa(x0[1],NamedTuple)
@test length(x0[1])==1
@test keys(x0[1])==(:x,)

@test_throws "basic component is undefined" rprocess(P,params=(a=1,k=3,x₀=5),x0=x0)
pomp!(P,rprocess=rlin)
x = rprocess(P,x0=x0,params=(a=1,k=3,x₀=5));
@test isa(x,Array{<:NamedTuple})
@test size(x)==(27,1,5)
@test length(x[17])==1
@test keys(x[17])==(:x,)

@test_throws "basic component is undefined" rmeasure(P,params=(a=1,k=3,x₀=5),x=x)
P = pomp(P,rmeasure=rmeas);
y = rmeasure(P,x=x,params=(a=1,k=3,x₀=5));
@test isa(y,Array{<:NamedTuple})
@test size(y)==(27,1,5)
@test length(y[37])==1
@test keys(y[37])==(:pop,)

@test_throws "basic component is undefined" dmeasure(P,x=x,params=(a=1,k=3,x₀=5))
pomp!(P,dmeasure=dmeas,params=(a=1,k=3,x₀=5))
ell = dmeasure(P,y=y[:,:,[1]],x=x);
@test isa(ell,Array{<:Real,3})
@test size(ell)==(27,1,5)
simulate!(P,dmeasure=dmeas,params=(a=1,k=3,x₀=5))
ell = dmeasure(P)
@test size(ell)==(27,1,1)

pomp!(P,rprocess=nothing);
@test_throws "basic component is undefined" rprocess(P,params=(a=1,k=3,x₀=5),x0=x0)

@test_throws "is not defined" POMP.val_array("yes")
@test size(POMP.val_array(y,3,3))==(3,3,15)
@test_throws "size mismatch" POMP.val_array(y,3,2)

pomp!(P,rinit=function(;_...) error("yikes!") end)
@test_throws "in `rinit`: yikes!" rinit(P,params=(x₀=3,))
pomp!(P,rmeasure=function(;_...) error("yikes!") end)
@test_throws "in `rmeasure`: yikes!" rmeasure(P,params=(a=1,),x=x)
pomp!(P,dmeasure=function(;_...) error("yikes!") end)
@test_throws "in `dmeasure`: yikes!" dmeasure(P,params=(a=1,),x=x)
pomp!(P,rprocess=function(;_...) error("yikes!") end)
@test_throws "in `rprocess`: yikes!" rprocess(P,params=(x₀=3,),x0=x0)

