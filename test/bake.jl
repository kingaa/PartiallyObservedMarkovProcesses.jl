using PartiallyObservedMarkovProcesses
import PartiallyObservedMarkovProcesses as POMP
using Random: seed!, default_rng, Xoshiro
using Test

@info h1("testing reproducibility functions")

@testset verbose=true "reproducibility" begin

    f(x) = x+11
    g(x) = 3*x*rand()

    tmpfile = tempname(suffix=".jlso")

    z1 = @bake tmpfile begin
        f(3)+g(2)
    end

    z2 = @bake tmpfile begin
        f(3)+g(2)
    end
    @test z1 == z2

    z3 = @test_logs (:warn,r"recomputing") begin
        @bake tmpfile begin
            f(3)+g(3)
        end
    end
    @test z1 != z3

    r = default_rng()
    seed!(r,58361444)
    x1 = rand(r,2)
    y1 = @freeze r 39 rand(r,2)
    x2 = rand(r,2)
    seed!(r,58361444)
    x3 = rand(r,4)
    y2 = @freeze r 39 rand(r,2)
    @test [x1...,x2...] == x3
    @test y1==y2
    seed!(r,58361444)
    y3 = @freeze 39 rand(r,2)
    @test y1==y3

    r = Xoshiro()
    seed!(r,58361444)
    x1 = rand(r,2)
    y1 = @freeze r 39 rand(r,2)
    x2 = rand(r,2)
    seed!(r,58361444)
    x3 = rand(r,4)
    y2 = @freeze r 39 rand(r,2)
    @test [x1...,x2...] == x3
    @test y1==y2

    r = default_rng()
    seed!(r,58361444)
    x1 = rand(r,2)
    z1 = @test_logs (:warn,r"recomputing") begin
        @bake tmpfile begin
            @freeze r 33 f(3)+g(2)
        end
    end
    x2 = rand(r,2)

    seed!(r,58361444)
    x3 = rand(r,2)
    z2 = begin
        @bake tmpfile begin
            @freeze r 33 f(3)+g(2)
        end
    end
    x4 = rand(r,2)

    seed!(r,58361444)
    x5 = rand(r,2)
    z3 = @test_logs (:warn,r"recomputing") begin
        @freeze r 33 @bake tmpfile begin
            f(3)+g(2)
        end
    end
    x6 = rand(r,2)

    seed!(r,58361444)
    z4 = begin
        @freeze r 33 @bake tmpfile begin
            f(3)+g(2)
        end
    end

    seed!(r,58361444)
    x7 = rand(r,2)
    x8 = rand(r,2)

    @test x1==x3
    @test x2==x4
    @test x1==x5
    @test x2==x6
    @test x1==x7
    @test x2==x8
    @test z1==z2
    @test z1==z3

end
