using POMP
using Random: seed!
using Test

@testset verbose=true "constructor and workhorses" begin

    seed!(263260083)

    include("test1.jl")

    @test_throws "must be no later than" pomp(times=0:20,t0=5)
    @test_throws "times must be nondecreasing" pomp(t0=0,times=[-t for t in 0:7])
    @test_throws "same elementary type" pomp(times=0:20,t0=0.0)

    P = pomp(times=0:20,t0=0);
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

    P = pomp(P,rprocess=discrete_time(rlin,dt=1));
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

end
