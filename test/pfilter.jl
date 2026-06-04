using PartiallyObservedMarkovProcesses
import PartiallyObservedMarkovProcesses as POMP
using Distributions
using Random
using Test
using BenchmarkTools

@info h1("pfilter tests")

@testset verbose=true "pfilter" begin

    Random.seed!(263260083)

    rin = function(;x0,_...)
        d = Poisson(x0)
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

    p1 = (a=1.5,k=7.0,x0=5.0);

    P = simulate(
        t0=0,
        times=0:20,
        params=p1,
        rinit=rin,
        rprocess=discrete_time(rlin,dt=1),
        rmeasure=rmeas,
        logdmeasure=logdmeas
    );
    P = P[1];
    @test_throws r"keyword argument .* not assigned" simulate(P,params=(a=1.5,k=7.0))

    P = pomp(
        obs(P),
        times=times(P),
        t0=0,
        rinit=rin,
        rprocess=discrete_time(rlin),
        rmeasure=rmeas,
        logdmeasure=logdmeas
    );
    @test P isa POMP.PompObject

    x0 = rinit(P,params=p1,nsim=10);
    y = obs(P);
    t = times(P);
    x = similar(x0,1,size(x0)...);
    rprocess!(P,x,x0=x0,times=t[1:1],params=[p1])
    @test x0==x[1,:,:]

    Q = pfilter(P,Np=1000,params=p1);
    @test :x ∈ paramsymbs(Q)
    @test Q isa POMP.PfilterdPompObject
    @test occursin(r"PfilterdPompObject .* Np=",sprint(show,Q))
    @test all(Q.x0.==Q.pred[1,:])
    @test_throws r"keyword argument .* not assigned" pfilter(Q,params=(a=1.5,k=7.0));
    @btime pfilter($Q,params=(a=1.5,k=7.0,x0=5.0));
    @btime pfilter($Q,params=(k=7.0,a=1.5,x0=5.0));
    x0 = rinit(Q,nsim=5)
    @test x0 isa Array{<:NamedTuple}
    @test size(x0)==(1,5)
    rinit!(Q,x0)

    d = melt(Q);
    @test size(d)==(21,4)
    @test propertynames(d)==[:time; :y; :ess; :cond_logLik]

    P1 = pomp(P;logdmeasure=function (;_...) -Inf end)
    Q1 = pfilter(P1,Np=100,params=p1)
    @test isinf(logLik(Q1))
    @test all(eff_sample_size(Q1).==0)
    @test all(eff_sample_size(Q1).==0)
    @test all(isinf.(cond_logLik(Q1)))
    @test Q1.filt==Q1.pred

    Q2 = [pfilter(Q,Np=100) for _ ∈ 1:5]
    @test logmeanexp(logLik.(Q2)) isa Float64
    @test logmeanexp(logLik.(Q2),se=true) isa @NamedTuple{est::Float64,se::Float64}
    @test logmeanexp(logLik.(Q2),ess=true) isa @NamedTuple{est::Float64,ess::Float64}
    @test logmeanexp(logLik.(Q2),ess=true,se=true) isa @NamedTuple{est::Float64,se::Float64,ess::Float64}
    @test logmeanexp(logLik.(Q2)) > mean(logLik.(Q2))

end
