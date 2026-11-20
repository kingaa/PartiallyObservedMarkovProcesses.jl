using PartiallyObservedMarkovProcesses
using Test
using DataFrames: DataFrame

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
    @test obs(P) isa Vector{@NamedTuple{}}
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
    @test 8.0 == @with coef(Q[1]) k+a

end

@testset verbose=true "@with" begin

    p = (β=16,γ=2,ω=0,ψ=0,χ=1,pop=100,S0=30,I0=10,R0=60)

    trans = @with p begin
        m = inv(S0+I0+R0)
        R0 = β/γ
        S = S0/pop
        (;R0,S)
    end
    @test trans==(R0=8,S=0.3,)

    x = DataFrame(:a=>1:5, :b=>23)
    @test @with(x,5) == 5
    y = @with x a+b
    @test y == 24:28
    g = @with x function(x) x.+(a.+b) end
    @test g(3) == 27:31
    @test_throws "will not work" @with y 5

end
