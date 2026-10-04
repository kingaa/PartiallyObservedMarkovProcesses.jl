using PartiallyObservedMarkovProcesses
import PartiallyObservedMarkovProcesses as POMP
using PartiallyObservedMarkovProcesses.Examples
using Distributions: LogNormal, logpdf
using Random
using DataFrames
using Test

@info h1("LogBaryNormal perturbations, and mif without a declared init_state")

# Two fixes made at 4f743d8 on 2026-09-22 are now upstream (027c7bc): the
# `barycentric` normalization and `mif` on a model without `init_state`.
# Upstream's tests cover `barycentric` (test/transformations.jl) and the
# average behaviour of `LogBaryNormal` (test/perturbn.jl); this file keeps
# what they do not assert: every draw on the simplex, exact projection at
# zero scale, the returned names, exact behaviour at the other lags, and
# `mif` on a model without a declared `init_state`.

@testset verbose=true "core fixes" begin

    @testset "LogBaryNormal keeps the group on the simplex" begin
        Random.seed!(3301)
        # one standard deviation for the whole group (upstream's reworked
        # @perturbn no longer offers one per member)
        f1 = @perturbn (g,h,i)~LogBaryNormal(3.0)
        for f ∈ (f1,), scale ∈ (0.01,1.0)
            for start ∈ ((g=0.2,h=0.3,i=0.5),(g=2.0,h=3.0,i=5.0))
                draws = [f(scale,1;start...) for _ ∈ 1:5000]
                @test all(isapprox(d.g+d.h+d.i,1.0) for d ∈ draws)
                @test all(min(d.g,d.h,d.i) ≥ 0 for d ∈ draws)
            end
        end
        # a group that starts off the simplex is projected onto it, even
        # with no perturbation at all
        z = f1(0.0,1;g=2.0,h=3.0,i=5.0)
        @test collect(z) ≈ [0.2,0.3,0.5]
        # the draw is random: two draws differ
        @test f1(1.0,1;g=0.2,h=0.3,i=0.5) != f1(1.0,1;g=0.2,h=0.3,i=0.5)
        # other parameters are untouched and not returned
        f3 = @perturbn (g,h)~LogBaryNormal(1.0)
        d = f3(1.0,1;g=0.4,h=0.6,other=7.0)
        @test keys(d) == (:g,:h) && d.g+d.h ≈ 1
        # as an initial-value perturbation: applied at lag 0 only
        f4 = @perturbn (g,h,i)~ivp(LogBaryNormal(1.0))
        d0 = f4(1.0,0;g=0.2,h=0.3,i=0.5)
        @test d0 != (g=0.2,h=0.3,i=0.5) && d0.g+d0.h+d0.i ≈ 1
        @test f4(1.0,3;g=0.2,h=0.3,i=0.5) == (g=0.2,h=0.3,i=0.5)
    end

    @testset "multi-iteration mif without a declared init_state" begin
        # The SIR example as shipped declares init_state; the copy below
        # declares none, while its rinit still returns a four-component
        # integer state (S,I,R,C), C being an accumulator. This is how
        # PhyloPOMP builds its genealogy filters.
        P = sir(times=range(start=1.0,stop=20.0,step=1.0))
        Q = pomp(P;init_state=nothing)
        @test init_state(Q) == (;)
        θ = merge(coef(P),(β=1.0,γ=0.2))
        x0 = rinit(Q;params=θ)[1]
        @test x0 isa @NamedTuple{S::Int64,I::Int64,R::Int64,C::Int64}
        @test x0.S+x0.I+x0.R == θ.N
        ptb = @perturbn(
            β ~ LogNormal(0.05),
            γ ~ LogNormal(0.05),
            (S0,I0,R0) ~ ivp(LogBaryNormal(0.1)),
        )
        Nmif = 6
        fit(M;kw...) = begin
            Random.seed!(4417)
            mif(M;params=θ,Np=200,Nmif,perturbations=ptb,
                cooling=geometric_cooling(0.5),avfun=x->sum(x)/length(x),kw...)
        end
        for (trigger,target) ∈ ((missing,missing),(0.5,0.5))
            A = fit(P;trigger,target)
            B = fit(Q;trigger,target)
            @test B isa POMP.MifdPompObject
            tB = traces(B)
            # the starting point and one row per iteration (numbered from 1,
            # as `traces` does), all likelihoods finite
            @test tB.iteration == 1:Nmif+1
            @test all(isfinite,skipmissing(tB.logLik))
            @test count(!ismissing,tB.logLik) == Nmif
            # the estimates move
            @test tB.β[end] != θ.β && tB.γ[end] != θ.γ
            # the final filter carries the states rinit returns
            @test eltype(B.pfobj.filt) == @NamedTuple{S::Int64,I::Int64,R::Int64,C::Int64}
            @test isfinite(logLik(B))
            # the arithmetic mean of points on the simplex is on the simplex
            @test all(isapprox.(tB.S0 .+ tB.I0 .+ tB.R0,1.0))
            # identical to the model that declares its state
            @test isequal(tB,traces(A))
            @test coef(B) == coef(A)
            @test logLik(B) == logLik(A)
            # continuation, several more iterations
            Random.seed!(4418); A2 = mif(A;Nmif=3)
            Random.seed!(4418); B2 = mif(B;Nmif=3)
            @test traces(B2).iteration == 1:4
            @test isequal(traces(B2),traces(A2))
        end
    end

end
