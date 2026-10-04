using PartiallyObservedMarkovProcesses
import PartiallyObservedMarkovProcesses as POMP
using DataFrames
using Distributions
using Random
using Test

@info h1("profile tests")

@testset verbose=true "profile" begin

    Random.seed!(1832445720)

    rin = function(;x0,_...)
        (x=rand(Poisson(x0)),)
    end
    rlin = function (;t,a,x,_...)
        (x=rand(Poisson(a*x)),)
    end
    rmeas = function (;x,k,_...)
        (y=rand(NegativeBinomial(k,k/(k+x))),)
    end
    logdmeas = function (;x,y,k,_...)
        logpdf(NegativeBinomial(k,k/(k+x)),y)
    end

    p1 = (a=1.5,k=7.0,x0=5.0)

    P = simulate(
        t0=0,
        times=0:20,
        params=p1,
        rinit=rin,
        rprocess=discrete_time(rlin,dt=1),
        rmeasure=rmeas,
        logdmeasure=logdmeas
    )[1]

    ptb = @perturbn(k ~ LogNormal(0.05), x0 ~ LogNormal(0.05))
    cool = geometric_cooling(0.5)
    d = profile_design(
        a=[1.2,1.35,1.5,1.65,1.8];
        lower=(k=3.0,x0=3.0),upper=(k=10.0,x0=8.0),
        nprof=2,rng=MersenneTwister(5),
    )

    @testset "profile" begin
        Random.seed!(33)
        pr = profile(P,d;Nmif=3,Np=50,perturbations=ptb,cooling=cool,nreps=2,Np_eval=100)
        @test pr isa DataFrame
        @test nrow(pr) == 10
        @test propertynames(pr) == [:a,:k,:x0,:loglik,:se,:ess]
        # the profiled parameter is untouched, the others have moved
        @test pr.a == d.a
        @test all(pr.k .!= d.k)
        @test all(pr.x0 .!= d.x0)
        @test all(pr.k .> 0) && all(pr.x0 .> 0)
        @test all(isfinite,pr.loglik)
        @test all(pr.ess .≤ 2)
        # reproducible under the same seed
        Random.seed!(33)
        pr2 = profile(P,d;Nmif=3,Np=50,perturbations=ptb,cooling=cool,nreps=2,Np_eval=100)
        @test pr2 == pr
        # a grid written with integers, e.g. `a = 1:2`, is a grid of real
        # values: `mif`'s default geometric mean need not return an integer
        # exactly, so an integer column would fail
        di = profile_design(a=1:2;lower=(k=3.0,x0=3.0),upper=(k=10.0,x0=8.0),nprof=1,rng=MersenneTwister(6))
        pri = profile(P,di;Nmif=1,Np=20,perturbations=ptb,cooling=cool)
        @test pri.a == [1.0,2.0] && all(isfinite,pri.loglik)
        # mcap runs on the profile output
        m = mcap(pr.loglik,pr.a;span=1.0)
        @test m isa MCAP
        @test minimum(pr.a) ≤ m.mle ≤ maximum(pr.a)
    end

    @testset "fixed parameters" begin
        # perturbing a profiled parameter is refused
        bad = @perturbn(a ~ LogNormal(0.1), k ~ LogNormal(0.1))
        @test_throws r"must not be perturbed" profile(P,d;Nmif=1,Np=10,perturbations=bad,cooling=cool)
        badivp = @perturbn(k ~ LogNormal(0.1), a ~ ivp(LogNormal(0.1)))
        @test_throws r"must not be perturbed" profile(P,d;Nmif=1,Np=10,perturbations=badivp,cooling=cool)
        # the record survives row subsetting of the design
        @test_throws r"must not be perturbed" profile(P,d[d.a .< 1.5,:];Nmif=1,Np=10,perturbations=bad,cooling=cool)
        # a plain data frame without metadata is accepted as a design,
        # with a warning that nothing is protected; an unperturbed
        # parameter keeps its exact value
        pk = @perturbn(k ~ LogNormal(0.05))
        pr3 = @test_logs (:warn,r"no parameter is protected") profile(P,DataFrame(a=[1.5],k=[7.0],x0=[5.0]);
            Nmif=1,Np=20,perturbations=pk,cooling=cool)
        @test nrow(pr3) == 1
        @test pr3.x0 == [5.0] && pr3.a == [1.5]
        # parameters absent from the design appear in the output
        pr4 = profile(P,DataFrame(a=[1.5]);Nmif=1,Np=20,perturbations=pk,cooling=cool,profiled=(:a,))
        @test propertynames(pr4) == [:a,:k,:x0,:loglik,:se,:ess]
        @test pr4.x0 == [5.0] && pr4.a == [1.5] && pr4.k[1] != 7.0
    end

    @testset "the profiled keyword" begin
        # `hcat` drops the record that `profile_design` keeps
        dh = hcat(d[:,[:a,:k]],DataFrame(x0=d.x0))
        @test !haskey(metadata(dh),"profiled")
        bad = @perturbn(a ~ LogNormal(0.1), k ~ LogNormal(0.1))
        pk = @perturbn(k ~ LogNormal(0.05))
        # naming the profiled parameter protects it again
        @test_throws r"must not be perturbed" profile(P,dh;Nmif=1,Np=10,perturbations=bad,cooling=cool,profiled=(:a,))
        @test_throws r"must not be perturbed" profile(P,dh;Nmif=1,Np=10,perturbations=bad,cooling=cool,profiled=:a)
        # without it, a single warning says that nothing is protected
        @test_logs (:warn,r"no parameter is protected") profile(P,dh[1:3,:];Nmif=1,Np=10,perturbations=pk,cooling=cool)
        # `profiled = ()` declares that nothing is profiled: no warning
        @test_logs profile(P,dh[1:2,:];Nmif=1,Np=10,perturbations=pk,cooling=cool,profiled=())
        # the record and the keyword together protect both
        @test_throws r"`k` must not be perturbed" profile(P,d;Nmif=1,Np=10,perturbations=pk,cooling=cool,profiled=(:k,))
        @test_throws r"`bogus` is not a parameter" profile(P,d;Nmif=1,Np=10,perturbations=pk,cooling=cool,profiled=(:bogus,))
    end

    @testset "params supplies parameters the design lacks; the design wins" begin
        # the perturbation function records the fixed value of `a` it is
        # given: the fit must see the design's values, not those of `params`
        pk = @perturbn(k ~ LogNormal(0.05))
        # (the rows may run on several threads, so the record is locked)
        seen = Float64[]; lk = ReentrantLock()
        rec = (scale, lag; k, a, _...) -> (lock(() -> push!(seen,a),lk); (k=k*exp(0.05*scale*randn()),))
        Random.seed!(47)
        pr = profile(P,DataFrame(a=[1.2,1.8]);Nmif=1,Np=10,perturbations=rec,cooling=cool,
            profiled=(:a,),params=(a=5.0,))
        @test pr.a == [1.2,1.8]
        @test 5.0 ∉ seen && Set(seen) == Set([1.2,1.8])
        # a model that stores no parameters gets them from `params`
        Q = pomp(P;params=nothing)
        pq = profile(Q,DataFrame(a=[1.5]);Nmif=1,Np=10,perturbations=pk,cooling=cool,
            profiled=(:a,),params=p1)
        @test pq.a == [1.5] && pq.x0 == [5.0] && isfinite(pq.loglik[1])
        # with `params` holding only what the design lacks
        pq2 = profile(Q,DataFrame(a=[1.5]);Nmif=1,Np=10,perturbations=pk,cooling=cool,
            profiled=(:a,),params=(k=7.0,x0=5.0))
        @test pq2.a == [1.5] && pq2.x0 == [5.0] && isfinite(pq2.loglik[1])
        # reserved names are refused in `params` too
        @test_throws r"`se` is reserved" profile(P,DataFrame(a=[1.5]);Nmif=1,Np=10,perturbations=pk,
            cooling=cool,profiled=(:a,),params=(se=1.0,))
    end

    @testset "model components passed to mif are used in the evaluation" begin
        # with a flat measurement density of -1 at each of 21 times, the
        # log likelihood is exactly -21 whatever the parameters
        flat = function (;_...) -1.0 end
        pr = profile(P,DataFrame(a=[1.5]);Nmif=1,Np=20,
            perturbations=@perturbn(k ~ LogNormal(0.05)),cooling=cool,
            logdmeasure=flat,profiled=(:a,))
        @test pr.loglik[1] ≈ -length(times(P))
    end

    @testset "evaluation uses the resampling settings of the fit" begin
        # the resampling settings are honored.  In the frozen-state model
        # with 2 particles, a likelihood estimate of 0.5 or 1.0 can only
        # arise if the filter resampled; with trigger = 0 it never does,
        # so every estimate is 0, 0.9, or 1.8
        Pz = pomp([(y=0.0,),(y=0.0,)]; t0=0.0, times=[1.0,2.0], params=(d=1.0,),
            rinit=(;_...) -> (x = rand() < 0.5 ? 1.0 : 0.0,),
            rprocess=discrete_time((;x,_...) -> (x=x,),dt=1.0),
            logdmeasure=(;t,x,_...) -> t == 1 ? log(x == 1 ? 1.8 : 0.2) : log(x == 1 ? 1.0 : 0.0))
        Random.seed!(46)
        pr = profile(Pz,DataFrame(d=ones(40));Nmif=2,Np=10,Np_eval=2,
            perturbations=@perturbn(d ~ LogNormal(0.1)),cooling=cool,
            avfun=x->sum(x)/length(x),trigger=0.0,target=0.0,profiled=())
        @test round.(exp.(pr.loglik),digits=6) ⊆ [0.0,0.9,1.8]
        # and so are those carried by a `pfilter` or `mif` result
        pf = pfilter(Pz;Np=2,trigger=0.0,target=0.0)
        pr2 = profile(pf,DataFrame(d=ones(40));Nmif=1,Np=10,Np_eval=2,
            perturbations=@perturbn(d ~ LogNormal(0.1)),cooling=cool,
            avfun=x->sum(x)/length(x),profiled=())
        @test round.(exp.(pr2.loglik),digits=6) ⊆ [0.0,0.9,1.8]
    end

    @testset "evaluation is at exactly the reported point" begin
        # `c` is held fixed at 0.1, but the geometric mean of 0.1 repeated
        # is not exactly 0.1: an evaluation at the swarm average would pay
        # the penalty, -1000 at each of the 3 times
        Pc = pomp([(y=0.0,) for _ ∈ 1:3]; t0=0.0, times=[1.0,2.0,3.0], params=(c=0.1,k=1.0),
            rinit=(;_...)->(x=0.0,), rprocess=discrete_time((;x,_...)->(x=x,),dt=1.0),
            logdmeasure=(;c,_...) -> isequal(c,0.1) ? -1.0 : -1000.0)
        @test geomean(fill(0.1,20)) != 0.1
        pr = profile(Pc,DataFrame(c=[0.1]);Nmif=1,Np=20,cooling=cool,
            perturbations=@perturbn(k ~ LogNormal(0.1)),profiled=(:c,))
        @test pr.c == [0.1] && pr.loglik == [-3.0]
    end

    @testset "reserved names" begin
        pk = @perturbn(k ~ LogNormal(0.05))
        @test_throws r"`slice` is reserved" profile(pomp(P;params=merge(p1,(slice=1.0,))),DataFrame(a=[1.5]);
            Nmif=1,Np=10,perturbations=pk,cooling=cool)
        @test_throws r"`se` is reserved" profile(pomp(P;params=merge(p1,(se=1.0,))),DataFrame(a=[1.5]);
            Nmif=1,Np=10,perturbations=pk,cooling=cool)
    end

    @testset "finding the perturbed parameters" begin
        N = length(times(P))
        @test Set(POMP.perturbed_names(ptb,[p1],N)) == Set([:k,:x0])
        @test Set(POMP.perturbed_names(@perturbn(x0 ~ ivp(LogNormal(0.1))),[p1],N)) == Set([:x0])
        # the random numbers are left alone
        Random.seed!(44); u = rand()
        Random.seed!(44); POMP.perturbed_names(ptb,[p1],N); @test rand() == u
        # every lag `mif` uses is checked, not only the first two
        late = function (scale, lag; k, a, _...)
            k = rand(LogNormal(log(k),0.05*scale))
            lag == 7 ? (;k, a=rand(LogNormal(log(a),0.05*scale))) : (;k)
        end
        @test Set(POMP.perturbed_names(late,[p1],N)) == Set([:k,:a])
        @test_throws r"must not be perturbed" profile(P,d;Nmif=1,Np=10,perturbations=late,cooling=cool)
        # a parameter perturbed depending on its value, which the lags
        # cannot reveal, is caught when `mif` makes the call
        Random.seed!(45)
        sneaky = function (scale, lag; k, a, _...)
            move_a = k > 7.2    # decided by the incoming value of k
            k = rand(LogNormal(log(k),0.2*scale))
            move_a ? (;k, a=rand(LogNormal(log(a),0.05*scale))) : (;k)
        end
        @test POMP.perturbed_names(sneaky,[p1],N) == (:k,)
        @test_throws r"moved parameter `a`" profile(P,DataFrame(a=[1.5],k=[7.0],x0=[5.0]);
            Nmif=2,Np=20,perturbations=sneaky,cooling=geometric_cooling(1.0))
        # a function that moves `a` only below full scale, and moves it
        # back before each iteration ends: invisible both to probing and
        # to the trace
        sly = function (scale, lag; k, _...)
            k = rand(LogNormal(log(k),0.05*scale))
            scale == 1.0 ? (;k) : lag == 5 ? (;k, a=2.0) : lag == 15 ? (;k, a=1.5) : (;k)
        end
        @test POMP.perturbed_names(sly,[p1],N) == (:k,)
        @test_throws r"moved parameter `a`" profile(P,DataFrame(a=[1.5],k=[7.0],x0=[5.0]);
            Nmif=2,Np=20,perturbations=sly,cooling=geometric_cooling(0.5))
    end

end
