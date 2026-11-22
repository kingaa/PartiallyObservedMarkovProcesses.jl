using PartiallyObservedMarkovProcesses
using PartiallyObservedMarkovProcesses.Examples
using Distributions: LogNormal, logpdf
using Random
using DataFrames
using Test

@info h1("models without a declared init_state")

@testset verbose=true "undeclared init_state" begin

    # The Gompertz model twice: as shipped (init_state declared) and
    # built without init_state, as PhyloPOMP's genealogy filters are.
    # `mif` used to fail on the second with a MethodError. Given the
    # same random numbers, the two must now give identical results.
    P = gompertz()
    Q = pomp(
        parus_data,t0=1960,times=:year,
        rinit=(;X0,_...) -> (;X=X0,),
        rprocess=discrete_time(
            function (;t,X,σₚ,r,K,_...)
                s = exp(-r)
                (;X=rand(LogNormal(s*log(X)+(1-s)*log(K),σₚ)),)
            end,
            dt=1,
        ),
        logdmeasure=(;pop,X,σₘ,_...) -> logpdf(LogNormal(log(X),σₘ),pop),
    )
    @test init_state(Q) == (;)
    p1 = (r=4.5,K=210.0,σₚ=0.7,σₘ=0.1,X0=150.0)
    ptb = @perturbn(K ~ LogNormal(0.1), σₚ ~ LogNormal(0.1))
    for (trigger,target) ∈ ((missing,missing),(0.5,0.5))
        fit(M) = begin
            Random.seed!(71)
            mif(M;params=p1,Np=100,Nmif=3,perturbations=ptb,
                cooling=geometric_cooling(0.5),trigger,target)
        end
        A = fit(P)
        B = fit(Q)
        @test isequal(traces(B),traces(A))
        @test coef(B) == coef(A)
        @test logLik(B) == logLik(A)
        # continuation
        Random.seed!(72); A2 = mif(A;Nmif=2)
        Random.seed!(72); B2 = mif(B;Nmif=2)
        @test isequal(traces(B2),traces(A2))
    end
    # and a profile runs on it
    pr = profile(Q,DataFrame(K=[200.0],r=[4.5],σₚ=[0.7],σₘ=[0.1],X0=[150.0]);
        Nmif=2,Np=50,perturbations=@perturbn(σₚ ~ LogNormal(0.1)),
        cooling=geometric_cooling(0.5))
    @test isfinite(pr.loglik[1]) && pr.K == [200.0]

end
