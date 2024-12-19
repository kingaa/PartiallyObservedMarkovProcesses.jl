using PartiallyObservedMarkovProcesses
using Test

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
    @test isa(times(P),Vector{<:Real})
    @test size(times(P))==(21,)
    @test isa(init_state(P),Nothing)
    @test isa(init_state(Q[1]),NamedTuple)
    @test isa(init_state(Q),Array{<:NamedTuple,2})
    @test size(init_state(Q))==(2,3)
    @test isa(obs(P),Nothing)
    @test isa(obs(Q[1]),Array{<:NamedTuple,1})
    @test size(obs(Q[1]))==(21,)
    @test isa(obs(Q),Array{<:NamedTuple,3})
    @test size(obs(Q))==(21,2,3)
    @test isa(states(P),Nothing)
    @test isa(states(Q[2]),Array{<:NamedTuple,1})
    @test size(states(Q[2]))==(21,)
    @test isa(states(Q),Array{<:NamedTuple,3})
    @test size(states(Q))==(21,2,3)
    @test isa(coef(P),Nothing)
    @test isa(coef(Q[3]),NamedTuple)
    @test isa(coef(Q),Array{<:NamedTuple,2})
    @test size(coef(Q))==(2,3)
    @test isa(coef(Q[3],:k,:a),NamedTuple)
    @test isa(coef(Q,:k,:a),Array{<:NamedTuple,2})
    
end
