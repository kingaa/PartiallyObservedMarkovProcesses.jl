using PartiallyObservedMarkovProcesses
import PartiallyObservedMarkovProcesses as POMP
using PartiallyObservedMarkovProcesses.Examples
using Random
using DataFrames
using Distributions: LogNormal
using Test
using BenchmarkTools

@info h1("mif tests")

@testset verbose=true "mif" begin

    Random.seed!(263260083)

    P = sir(times=range(start=1.0,stop=10.0,step=1.0))
    Pf = pfilter(P,Np=1000,trigger=0.2,target=0.5)

    pkern(s, lag; β, γ, S0, I0, R0, _...) = begin
        β = rand(LogNormal(log(β),0.05*s))
        γ = rand(LogNormal(log(γ),0.05*s))
        if lag == 0
            S0 = rand(LogNormal(log(S0),0.3*s))
            I0 = rand(LogNormal(log(I0),0.3*s))
            R0 = rand(LogNormal(log(R0),0.3*s))
            m = S0 + I0 + R0
            S0 /= m
            I0 /= m
            R0 /= m
        end
        (;β,γ,S0,I0,R0,)
    end

    cool(k) = begin
        0.5^(Float64(k)/50)
    end

    p1 = merge(coef(P),(β=1.0,γ=0.2,S0=0.8))
    Pf0 = pfilter(Pf,params=p1)
    
    M = mif(
        P,params=p1,
        Np=1000,Nmif=200,
        perturbation_kernel=pkern,
        cooling_schedule=cool,
    )
    @test M isa POMP.MifdPompObject
    M2 = mif(M,Nmif=2)
    @test M2 isa POMP.MifdPompObject
    @test M2.pfobj.Np == M.pfobj.Np

    @test traces(M) isa DataFrame
    @test all(propertynames(coef(M)) .∈ Ref(propertynames(traces(M))))
    @test :logLik ∈ propertynames(traces(M))
    @test :iteration ∈ propertynames(traces(M))

    @test logLik(M) > logLik(Pf) > logLik(Pf0)

    @test melt(coef(M)) == DataFrame(traces(M)[end,[keys(coef(M))...]])
    @test pomp(M) isa POMP.PompObject
    @test pfilter(M) isa POMP.PfilterdPompObject
    @test_throws "Incorrect call" mif("hello!")

end
