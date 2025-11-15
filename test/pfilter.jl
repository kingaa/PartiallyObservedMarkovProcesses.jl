using PartiallyObservedMarkovProcesses
import PartiallyObservedMarkovProcesses as POMP
using Distributions
using Random
using Test

println("- pfilter tests")

@testset verbose=true "pfilter" begin

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

    p1 = (a=1.5,k=7.0,x₀=5.0);

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
    @test_throws "in `simulate`" simulate(P,params=(a=1.5,k=7.0))

    P = pomp(
        obs(P),
        times=times(P),
        t0=0,
        rinit=rin,
        rprocess=discrete_time(rlin),
        rmeasure=rmeas,
        logdmeasure=logdmeas
    );
    @test isa(P,POMP.PompObject)

    x0 = rinit(P,params=p1,nsim=10);
    y = obs(P);
    t = times(P);
    x = similar(x0,1,size(x0)...);
    rprocess!(P,x,x0=x0,times=t[1:1],params=[p1])
    @test x0==x[1,:,:]

    Q = pfilter(P,Np=1000,params=p1);
    @test isa(Q,POMP.PfilterdPompObject)
    println("    ",Q)
    @time pfilter(Q,params=(a=1.5,k=7.0,x₀=5.0));
    @time pfilter(Q,params=(a=1.5,k=7.0,x₀=5.0));
    @time pfilter(Q,params=(k=7.0,a=1.5,x₀=5.0));
    @time pfilter(Q,params=(k=7.0,a=1.5,x₀=5.0));
    @test all(Q.x0.==Q.pred[1,:])
    @test_throws "in `pfilter`: in `rinit`" pfilter(Q,params=(a=1.5,k=7.0));

    d = melt(Q);
    @test size(d)==(21,4)
    @test propertynames(d)==[:time; :y; :ess; :cond_logLik]

    P = pomp(P;logdmeasure=function (;_...) -Inf end)
    Q = pfilter(P,Np=100,params=p1)
    @test isinf(Q.logLik)
    @test all(Q.eff_sample_size.==0)
    @test all(Q.eff_sample_size.==0)
    @test all(isinf.(Q.cond_logLik))
    @test Q.filt==Q.pred

end
