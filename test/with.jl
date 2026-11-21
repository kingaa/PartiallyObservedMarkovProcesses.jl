using PartiallyObservedMarkovProcesses
using Test
using DataFrames: DataFrame

@info h1("testing @with")

@testset verbose=true "@with" begin

    S = 100
    R0 = 0
    p = (β=16,γ=2,ω=0,ψ=0,χ=1,pop=100,S0=30,I0=10,R0=60)
    trans = @with p begin
        m = inv(S0+I0+R0)
        R0 = β/γ
        S = S0/pop
        (;R0,S)
    end
    @test trans==(R0=8,S=0.3,)
    @test R0 == 0 && S == 100

    x = DataFrame(:a=>1:5, :b=>23)
    @test @with(x,5) == 5

    y = @with x a+b
    @test y == 24:28

    g = @with x function(x) x.+(a.+b) end
    @test g(3) == 27:31
    @test_throws "will not work" @with y 5

end
