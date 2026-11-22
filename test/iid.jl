using PartiallyObservedMarkovProcesses
import PartiallyObservedMarkovProcesses as POMP
using Distributions
using Random
using Statistics
using Test

@info h1("pfilter tests on iid model")

@testset verbose=true "iid" begin

    Random.seed!(263260083)

    rin = function(;σ,_...)
        d = Normal(0,1)
        (x=rand(d),)
    end

    rproc = function (;t,x,σ,_...)
        d = Normal(0,1)
        (x=rand(d),)
    end

    rmeas = function (;x,τ,_...)
        d = Normal(x,1)
        (y=rand(d),)
    end

    logdmeas = function (;x,y,τ,_...)
        d = Normal(x,τ)
        logpdf(d,y)
    end

    p1 = (σ=1.0,τ=1.0);

    P = simulate(
        t0=0,
        times=0:100,
        params=p1,
        rinit=rin,
        rprocess=discrete_time(rproc,dt=1),
        rmeasure=rmeas,
        logdmeasure=logdmeas
    )[1];

    loglik_iid(;y,σ,τ,_...) = begin
        d = Normal(0,sqrt(σ^2+τ^2))
        logpdf(d,y)
    end

    llexact = map(obs(P)) do x
        loglik_iid(;x...,p1...)
    end |> sum

    ll = logmeanexp([logLik(pfilter(P,Np=10000,params=p1,trigger=0.3,target=0.1)) for _ in 1:20],se=true,ess=true)
    @test abs(ll.est-llexact) < 4*ll.se

    # This latent process has no memory, so the test cannot detect a
    # missing normalizing constant at `target > 0`; the next one can.

end

@info h1("unbiasedness of the power-renormalized filter")

@testset verbose=true "power resampling is unbiased" begin

    # Frozen binary state X ~ Bernoulli(1/2), with g₁ = (0.2,1.8) and
    # g₂ = (0,1) at X = (0,1): the exact likelihood is 0.9.  With Np = 2
    # and target = 1/2, omitting the normalizing constant of the
    # tempered resampling step gives E[Ẑ] = 0.8875.

    Random.seed!(263260083)

    rin   = function(;_...); (x = rand() < 0.5 ? 1.0 : 0.0,); end
    rproc = function(;x,_...); (x=x,); end
    ldm   = function(;t,x,y,_...)
        t == 1 ? log(x == 1.0 ? 1.8 : 0.2) : log(x == 1.0 ? 1.0 : 0.0)
    end

    P = pomp(
        [(y=0.0,),(y=0.0,)];
        t0 = 0.0,
        times = [1.0,2.0],
        params = (dummy=0.0,),
        rinit = rin,
        rprocess = discrete_time(rproc,dt=1.0),
        logdmeasure = ldm,
    )

    zexact = 0.9
    nrep = 200_000

    # E[Ẑ] and its standard error, on the natural scale
    meanZ(;Np,trigger,target) = begin
        s = 0.0; ss = 0.0
        for _ ∈ 1:nrep
            z = exp(logLik(pfilter(P;Np,trigger,target)))
            s += z; ss += z*z
        end
        mu = s/nrep
        (est=mu, se=sqrt(max(ss/nrep-mu^2,0.0)/nrep))
    end

    # control: ordinary resampling, target = 0, must be unbiased
    r0 = meanZ(Np=2,trigger=1.0,target=0.0)
    @test abs(r0.est-zexact) < 4*r0.se

    for target ∈ [0.25,0.5,0.75]
        r = meanZ(Np=2,trigger=1.0,target=target)
        @info "target=$target: E[Ẑ]=$(round(r.est,digits=6)) ± $(round(r.se,digits=6)), exact=$zexact"
        @test abs(r.est-zexact) < 4*r.se
    end

    # reported, not tested: a bias would be O(1/Np)
    for Np ∈ [8,32,128]
        r = meanZ(Np=Np,trigger=1.0,target=0.5)
        @info "Np=$Np: E[Ẑ]=$(round(r.est,digits=6)) ± $(round(r.se,digits=6)), exact=$zexact"
    end

end
