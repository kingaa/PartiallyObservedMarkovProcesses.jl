using PartiallyObservedMarkovProcesses
using Test

@info h1("testing helper functions")

@testset verbose=true "helpers" begin

    include("test1.jl")

    P = pomp(
        times=1:21,
        t0=0,
        rinit=rin,
        rprocess=discrete_time(rlin),
        rmeasure=rmeas,
        logdmeasure=logdmeas
    )
    Q = simulate(P,params=[theta;theta],nsim=3);
    @test times(P) isa Vector{<:Real}
    @test size(times(P))==(21,)
    @test init_state(P) isa NamedTuple
    @test init_state(Q[1]) isa NamedTuple
    @test init_state(Q) isa Array{<:NamedTuple,2}
    @test size(init_state(Q))==(2,3)
    @test obs(P) isa Nothing
    @test obs(Q[1]) isa Array{<:NamedTuple,1}
    @test size(obs(Q[1]))==(21,)
    @test obs(Q) isa Array{<:NamedTuple,3}
    @test size(obs(Q))==(21,2,3)
    @test states(P) isa Nothing
    @test states(Q[2]) isa Array{<:NamedTuple,1}
    @test size(states(Q[2]))==(21,)
    @test states(Q) isa Array{<:NamedTuple,3}
    @test size(states(Q))==(21,2,3)
    @test isempty(coef(P))
    @test coef(Q[3]) isa NamedTuple
    @test coef(Q) isa Array{<:NamedTuple,2}
    @test size(coef(Q))==(2,3)
    @test coef(Q[3],:k,:a) isa NamedTuple
    @test coef(Q,:k,:a) isa Array{<:NamedTuple,2}

end
