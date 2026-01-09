using PartiallyObservedMarkovProcesses
import PartiallyObservedMarkovProcesses as POMP
using Test
using Random: seed!

import DifferentialEquations: AutoTsit5, Rosenbrock23

println("- vectorfield integration tests")

@testset "vectorfield integration" begin

    seed!(875002133)

    P = simulate(
        params=(
            ;r = 1, K = 1e4, A = 1e3,
            b = 1e-3, c = 1, m = 0.8,
            N₀ = 3000, P₀ = 4, t₀ = 0.0,
        ),
        t0=0.0,
        times=Float64.(0:10),
        rinit = function (;N₀,P₀,_...)
            (
                X=log(N₀),
                Y=log(P₀),
            )
        end,
        rmeasure = nothing,
        rprocess = vectorfield(
            function (X,Y;t,b,c,r,K,A,m,_...)
                N = exp(X)
                P = exp(Y)
                mu = [r*N, r*N*N/K, c*P*N/(1+N/A), m*P]
                dN = mu[1]-mu[2]-mu[3]
                dP = b*mu[3]-mu[4]
                [dN/N, dP/P]
            end,
            AutoTsit5(Rosenbrock23()),
            verbose=true,
            force_dtmin=true,
            dtmax=1.0
        ),
    )[1];

    @test isa(P,POMP.PompObject)
    @test isa(states(P),Array{@NamedTuple{X::Float64,Y::Float64},1})
    @test size(states(P))==(11,)

end
