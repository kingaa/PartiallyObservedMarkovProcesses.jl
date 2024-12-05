using POMP
using DataFrames
using Distributions
using Random
using RCall
using Test

@testset "basic construction and workhorses" begin

    Random.seed!(263260083)

    rin = function(;x₀,_...)
        d = Poisson(x₀)
        (x=rand(d),)
    end

    rlin = function (;t,a,x,_...)
        d = Poisson(a*x)
        (x=rand(d),)
    end

    rmeas = function (;x,k,_...)
        d = NegativeBinomial(k,k/(k+x))
        (y=rand(d),)
    end

    logdmeas = function (;x,y,k,_...)
        d = NegativeBinomial(k,k/(k+x))
        logpdf(d,y)
    end

    @test_throws "must be no later than" pomp(times=[t for t in 0:20],t0=5)
    @test_throws "times must be nondecreasing" pomp(t0=0,times=[-t for t in 0:7])
    @test_throws "same elementary type" pomp(times=[t for t in 0:20],t0=0.0)

    P = pomp(times=[t for t in 0:20],t0=0);
    @test isa(P,POMP.PompObject)
    @test timezero(P)==0
    @test isa(times(P),Vector{<:Real})
    @test length(times(P))==21
    @test isnothing(obs(P))
    @test isnothing(states(P))

    p1 = (a=1.0,k=7.0,x₀=5.0);
    p2 = (a=1.1,k=2.0,x₀=3.0);

    x0 = rinit(P,params=p1,nsim=3);
    @test isa(x0,Array{@NamedTuple{},2})
    @test size(x0)==(3,1)
    @test x0[1]==(;)
    rinit!(P,x0,params=p1);
    @test isa(x0,Array{@NamedTuple{},2})
    @test size(x0)==(3,1)
    @test x0[1]==(;)

    P = pomp(P,rinit=rin);
    x0 = rinit(P,params=[p1;p2],nsim=7);
    @test isa(x0,Array{<:NamedTuple,2})
    @test size(x0)==(7,2)
    @test keys(x0[8])==(:x,)
    @test_throws r"parameter .* undefined" rinit(P,params=(a=1.0,k=7.0))
    @test_throws r"parameter .* undefined" rinit!(P,x0,params=(a=1.0,k=7.0))

    x = rprocess(P,x0=x0,params=[p1;p2]);
    @test isa(x,Array{<:NamedTuple,3})
    @test size(x)==(7,2,21)
    @test x[:,:,21]==x0

    P = pomp(P,rprocess=euler(rlin,dt=1));
    x = rprocess(P,x0=x0,params=[p1;p2]);
    @test isa(x,Array{<:NamedTuple,3})
    @test size(x)==(7,2,21)
    @test x[:,:,1]==x0
    @test x[:,:,21]!=x0
    @test keys(x[17])==(:x,)
    @test_throws r"parameter .* undefined" rprocess(P,x0=x0[1,:],params=(k=7.0,x₀=5.0))

    y = rmeasure(P,x=x,params=[p1;p2]);
    @test isa(y,Array{@NamedTuple{},3})
    @test size(y)==(7,2,21)
    @test y[3]==(;)

    P = pomp(P,rmeasure=rmeas);
    y = rmeasure(P,x=x,params=[p1;p2]);
    @test isa(y,Array{<:NamedTuple,3})
    @test size(y)==(7,2,21)
    @test isassigned(y,3)
    @test keys(y[17])==(:y,)
    @test_throws r"parameter .* undefined" rmeasure(P,x=x[1,:,:],params=(a=1.0,))

    ell = logdmeasure(P,x=x,y=y,params=[p1;p2]);
    @test isa(ell,Array{POMP.LogLik,4})
    @test size(ell)==(size(y,1),size(x)...)
    @test all(ell.==0.0)

    P = pomp(P,logdmeasure=logdmeas);
    ell = logdmeasure(P,x=x,y=y,params=[p1;p2]);
    @test isa(ell,Array{POMP.LogLik,4})
    @test size(ell)==(size(y,1),size(x)...)
    @test all(ell.<=0)
    @test_throws r"parameter .* undefined" logdmeasure(P,y=y,x=x,params=(a=1.0,))

    @testset "val_array" begin
        @test POMP.val_array("yes")==["yes"]
        @test_throws "size mismatch" POMP.val_array(y,11,2)
        @test size(POMP.val_array(3,1,1,1))==(1,1,1)
        @test size(POMP.val_array(rand(3,2,2)))==(12,)
    end

    @testset "melt" begin
        d = melt(nothing);
        @test isa(d,DataFrame)
        @test size(d)==(0,0)
        d = melt((x=3,y=99));
        @test isa(d,DataFrame)
        @test size(d)==(1,2)
        @test names(d)==["x","y"]
        d = melt(fill((a=10,b=0,c=33),5,2));
        @test isa(d,DataFrame)
        @test size(d)==(10,5)
        @test names(d)==["id1","id2","a","b","c"]
        d = melt(fill((a=10,b=0,c=33),5,2),:i,:j);
        @test isa(d,DataFrame)
        @test size(d)==(10,5)
        @test names(d)==["i","j","a","b","c"]
    end

    @testset "simulate" begin

        p = [p1;p2];
        Q1 = simulate(P,params=p,nsim=3);
        d = leftjoin(
            melt(Q1),
            melt(p,:id2),
            on=:id2
        );
        R"""
library(tidyverse)
$d |>
  ggplot(aes(x=time,group=interaction(id1,id2),
    color=interaction(id1,id2)))+
  geom_line(aes(y=x))+
  geom_point(aes(y=y))+
  geom_line(aes(y=y),linetype=2)+
  guides(color="none")+
  facet_wrap(~a+k,ncol=1,scales="free_y",
    labeller=label_bquote(list(a==.(a),k==.(k))))+
  scale_y_continuous(transform="log1p")+
  theme_bw()
"""

        R"""ggsave(filename="test1-01.png",width=7,height=4)"""

    end

    @testset "helpers" begin
        Q = simulate(P,params=[p1;p2],nsim=3);
        @test isa(times(P),Vector{<:Real})
        @test size(times(P))==(21,)
        @test isa(init_state(P),Nothing)
        @test isa(init_state(Q[1]),NamedTuple)
        @test isa(init_state(Q),Array{<:NamedTuple,2})
        @test size(init_state(Q))==(3,2)
        @test isa(obs(P),Nothing)
        @test isa(obs(Q[1]),Array{<:NamedTuple,1})
        @test size(obs(Q[1]))==(21,)
        @test isa(obs(Q),Array{<:NamedTuple,3})
        @test size(obs(Q))==(21,3,2)
        @test isa(states(P),Nothing)
        @test isa(states(Q[2]),Array{<:NamedTuple,1})
        @test size(states(Q[2]))==(21,)
        @test isa(states(Q),Array{<:NamedTuple,3})
        @test size(states(Q))==(21,3,2)
        @test isa(coef(P),Nothing)
        @test isa(coef(Q[3]),NamedTuple)
        @test isa(coef(Q),Array{<:NamedTuple,2})
        @test size(coef(Q))==(3,2)
    end

    @testset "error messages" begin

        @test_throws "must be no later than" pomp(parus_data,times=:year,t0=1999)
        @test_throws "times must be nondecreasing" pomp(sort(parus_data,:pop),times=:year,t0=1940)

        P = pomp(P,rinit=function(;_...) error("yikes!") end);
        @test_throws "in `rinit`: yikes!" rinit(P,params=(x₀=3,))
        @test_throws "in `rinit!`: yikes!" rinit!(P,x0,params=(x₀=3,))
        P = pomp(P,rprocess=euler(function(;_...) error("yikes!") end,dt=1));
        @test_throws "in `rprocess!`: yikes!" rprocess(P,params=(x₀=3,),x0=x0)
        P = pomp(P,rmeasure=function(;_...) error("yikes!") end);
        @test_throws "in `rmeasure`: yikes!" rmeasure(P,params=(a=1,),x=x)
        P = pomp(P,logdmeasure=function(;_...) error("yikes!") end);
        @test_throws "in `logdmeasure!`: yikes!" logdmeasure(P,params=(a=1,),x=x,y=y)

        @test_throws "must be of the same length" pomp([(a=1,);(a=2,);(c=3,)];times=[1;2],t0=0)

    end
end
