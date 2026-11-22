using PartiallyObservedMarkovProcesses
using Random
using Test

@info h1("perturbations stay in their parameter spaces")

@testset verbose=true "perturbation support" begin

    # Each kernel must map its parameter space into itself, checked over
    # many draws at a scale large enough to leave the space were it not
    # preserved by construction.
    Random.seed!(88)
    # (LogBaryNormal and barycentric are tested in core_fixes.jl.)
    ex = @perturbn a~LogNormal(3.0) p~LogitNormal(3.0)
    draws = [ex(1.0,1;a=2.0,p=0.3) for _ ∈ 1:10000]
    @test all(d.a > 0 for d ∈ draws)
    @test all(0 < d.p < 1 for d ∈ draws)
    # a zero scale leaves every value where it was
    z = ex(0.0,1;a=2.0,p=0.3)
    @test z.a == 2.0 && z.p ≈ 0.3

end
