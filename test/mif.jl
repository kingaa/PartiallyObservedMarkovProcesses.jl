using PartiallyObservedMarkovProcesses
import PartiallyObservedMarkovProcesses as POMP
using PartiallyObservedMarkovProcesses.Examples
using Random
using DataFrames
using Distributions: LogNormal
using Test

@info h1("mif tests")

@testset verbose=true "mif" begin

    Random.seed!(263260083)

    P = sir(times=range(start=1.0,stop=10.0,step=1.0))
    Pf = pfilter(P,Np=1000,trigger=0.2,target=0.5)

    p1 = merge(coef(P),(β=1.0,γ=0.2,S0=0.8))
    Pf0 = pfilter(Pf,params=p1)

    M = mif(
        P,params=p1,Np=10,Nmif=2,
        perturbations=@perturbn(
            @lognormal(β,0.05),
            @lognormal(γ,0.05),
            @ivp(@logbarynormal((S0,I0,R0),0.3)),
        ),
        cooling=geometric_cooling(0.5)
    )
    @time M = mif(M,Np=1000,Nmif=200,)
    @test M isa POMP.MifdPompObject
    mif(M,Nmif=1,trigger=0.3,target=0.5)
    M2 = mif(M,Nmif=10,trigger=0.3,target=0.5)
    @test M2 isa POMP.MifdPompObject
    @test M2.pfobj.Np == M.pfobj.Np
    Pf1 = pfilter(M2,target=0.8)
    @test Pf1 isa POMP.PfilterdPompObject
    M2 = mif(
        Pf1,Nmif=10,trigger=0.1,
        perturbations=M.perturbations,
        cooling=geometric_cooling(0.1)
    )
    @test M2 isa POMP.MifdPompObject
    @test M2.pfobj.trigger==0.1
    @test M2.pfobj.target==0.8

    @test traces(M) isa DataFrame
    @test all(propertynames(coef(M)) .∈ Ref(propertynames(traces(M))))
    @test :logLik ∈ propertynames(traces(M))
    @test :iteration ∈ propertynames(traces(M))

    @test logLik(M) > logLik(Pf) > logLik(Pf0)

    @test melt(coef(M)) == DataFrame(traces(M)[end,[keys(coef(M))...]])
    @test pomp(M) isa POMP.PompObject
    @test pfilter(M) isa POMP.PfilterdPompObject
    @test eff_sample_size(M) isa Vector
    @test all(eff_sample_size(M) .≤ 1000)
    @test cond_logLik(M) isa Vector
    @test all(cond_logLik(M) .< 0)
    @test_throws "Incorrect call" mif("hello!")
    @test occursin(r"MifdPompObject .* Nmif=",sprint(show,M))

    ex = @perturbn @lognormal(a,0.1) @ivp(@lognormal(b,1),0) @logitnormal(p,0.1) @normal(c,10) @ivp(@normal(d,10),1) @ivp(@logbarynormal((e,f),1)) @logbarynormal((g,h,i),(1,2,3))

    x = ex(1, 0, a=1, b=10, c=0, d=3, p=0.8, e=1, f=3, g=3, h=1, i=1)
    @test x.b != 10 && x.d == 3
    x = ex(1, 1, p=0.8, a=1, b=10, c=0, d=7, e=1, f=3, g=3, h=1, i=1)
    @test x.b == 10 && x.d != 7 && x.e == 1
    x = ex(1, 10, a=1, b=10, c=0, p=0.8, d=5, e=1, f=3, g=3, h=1, i=1)
    @test x.b == 10 && x.d == 5

end
