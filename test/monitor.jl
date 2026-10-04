using PartiallyObservedMarkovProcesses
import PartiallyObservedMarkovProcesses as POMP
using PartiallyObservedMarkovProcesses.Examples
using DataFrames
using Random
using Test

@info h1("hyperbolic cooling, monitor, and resampled")

@testset verbose=true "mif additions" begin

    @testset "hyperbolic_cooling" begin
        c = hyperbolic_cooling(0.5)
        @test c(0) == 1
        @test c(50) ≈ 0.5
        a = 50*0.1/0.9
        @test hyperbolic_cooling(0.1)(7) ≈ a/(a+7)
        # slower than geometric: same end point at 50, more left later
        @test hyperbolic_cooling(0.1)(200) > geometric_cooling(0.1)(200)
        # `start` shifts the schedule for a continued computation
        c4 = hyperbolic_cooling(0.5;start=4)
        @test [c4(i) for i ∈ 0:5] == [c(i+4) for i ∈ 0:5]
        @test all(hyperbolic_cooling(1.0)(i) == 1 for i ∈ 0:10)
        @test_throws r"frac" hyperbolic_cooling(0.0)
        @test_throws r"frac" hyperbolic_cooling(1.5)
        @test_throws r"start" hyperbolic_cooling(0.5;start=-1)
    end

    P = gompertz()
    p1 = (r=4.5,K=210.0,σₚ=0.7,σₘ=0.1,X0=150.0)
    ptb = @perturbn(K ~ LogNormal(0.1), σₚ ~ LogNormal(0.1))
    Random.seed!(3)
    A = mif(P;params=p1,Np=200,Nmif=4,perturbations=ptb,cooling=hyperbolic_cooling(0.5))
    B = mif(A;Nmif=3,cooling=hyperbolic_cooling(0.5;start=4))

    @testset "resampled" begin
        @test resampled(A) isa Vector{Bool}
        @test length(resampled(A)) == length(times(A))
        @test resampled(A) == resampled(A.pfobj)
    end

    @testset "monitor" begin
        Random.seed!(9); u = rand(); Random.seed!(9)
        m = monitor([A,B];Np=300,nreps=2,every=2,seed=1)
        # the default random-number stream is left where it was
        @test rand() == u
        # numbered as in `traces` and counted on across the continuation:
        # the start, every 2nd iteration after it, and the last
        @test m.iteration == [1,3,5,7,8]
        @test propertynames(m) == [:iteration,:loglik,:se,:ess,:r,:K,:σₚ,:σₘ,:X0]
        @test all(isfinite,m.loglik)
        # the evaluation points are the trace's point estimates; the
        # continuation's first row, its perturbed starting point, is
        # omitted: it is not a completed iteration
        tA = traces(A); tB = traces(B)
        @test m.K[1:3] == tA.K[[1,3,5]]
        @test m.K[4:5] == tB.K[[3,4]]
        # reproducible from its seed, and a different seed differs
        @test monitor([A,B];Np=300,nreps=2,every=2,seed=1) == m
        @test monitor([A,B];Np=300,nreps=2,every=2,seed=2).loglik != m.loglik
        # each value is a fixed-parameter evaluation at its point
        Random.seed!(1)
        @test m.loglik == [
            pfilter_loglik(P;Np=300,nreps=2,params=NamedTuple(r[[:r,:K,:σₚ,:σₘ,:X0]])).loglik
            for r ∈ eachrow(m)
        ]
        # settings are recorded
        @test metadata(m,"Np") == 300 && metadata(m,"nreps") == 2
        @test metadata(m,"every") == 2 && metadata(m,"seed") == 1
        # a single run
        @test monitor(A;Np=100,seed=1).iteration == traces(A).iteration == 1:5
        @test_throws r"every" monitor(A;Np=100,seed=1,every=0)
        @test_throws r"no `mif` results" monitor(POMP.MifdPompObject[];Np=100,seed=1)
        # a parameter named `iteration` would overwrite the iteration column
        Mi = mif(P;params=merge(p1,(iteration=0.25,)),Np=20,Nmif=1,perturbations=ptb,cooling=hyperbolic_cooling(0.5))
        @test_throws r"`iteration` is reserved" monitor(Mi;Np=20,seed=1)
        # each point is evaluated under the model of the run that recorded
        # it: log densities -1 then -2 at each of 3 times give exactly -3, -6
        Q = pomp([(y=0.0,) for _ ∈ 1:3]; t0=0.0, times=[1.0,2.0,3.0], params=(a=1.0,),
            rinit=(;_...)->(x=0.0,), rprocess=discrete_time((;x,_...)->(x=x,),dt=1.0),
            logdmeasure=(;_...)->-1.0)
        C = mif(Q;Np=10,Nmif=2,perturbations=@perturbn(a ~ LogNormal(0.01)),cooling=geometric_cooling(0.5))
        D = mif(C;Nmif=2,logdmeasure=(;_...)->-2.0)
        @test monitor([C,D];Np=10,seed=1).loglik == [-3.0,-3.0,-3.0,-6.0,-6.0]
        # runs whose parameters differ are refused: a continuation that adds
        # `q` cannot be evaluated with the first run's parameters
        E = mif(C;Nmif=2,params=merge(coef(C),(q=0.1,)),logdmeasure=(;q=0.9,_...)->log(q))
        @test_throws r"differs in `q`" monitor([C,E];Np=10,seed=1)
        # and with that run's resampling settings.  In the frozen-state
        # model with 2 particles, a likelihood estimate of 0.5 or 1.0 can only
        # arise if the filter resampled; with trigger = 0 it never does,
        # so every estimate is 0, 0.9, or 1.8
        Pz = pomp([(y=0.0,),(y=0.0,)]; t0=0.0, times=[1.0,2.0], params=(d=1.0,),
            rinit=(;_...) -> (x = rand() < 0.5 ? 1.0 : 0.0,),
            rprocess=discrete_time((;x,_...) -> (x=x,),dt=1.0),
            logdmeasure=(;t,x,_...) -> t == 1 ? log(x == 1 ? 1.8 : 0.2) : log(x == 1 ? 1.0 : 0.0))
        Mz = mif(Pz;Np=10,Nmif=30,perturbations=@perturbn(d ~ LogNormal(0.1)),
            cooling=geometric_cooling(0.5),avfun=x->sum(x)/length(x),trigger=0.0,target=0.0)
        @test round.(exp.(monitor(Mz;Np=2,seed=1).loglik),digits=6) ⊆ [0.0,0.9,1.8]
    end

end
