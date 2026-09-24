using PartiallyObservedMarkovProcesses
import PartiallyObservedMarkovProcesses as POMP
using Test

@info h1("parameter transformations")

@testset verbose=true "parameter transformations" begin

    @testset "barycentric on tuples and NamedTuples" begin
        ## a tuple is projected onto the unit simplex, keeping its type
        x = (2.0,3.0,5.0)
        b = barycentric(x)
        @test b isa NTuple{3,Float64}
        @test collect(b) ≈ [0.2,0.3,0.5]
        @test sum(b) ≈ 1
        ## onto the simplex of total n
        @test collect(barycentric(x,10)) ≈ [2.0,3.0,5.0]
        @test collect(barycentric((0.9,0.01,0.1),10000)) ≈ 10000 .* [0.9,0.01,0.1] ./ 1.01
        @test sum(barycentric((0.9,0.01,0.1),10000)) ≈ 10000
        ## a NamedTuple keeps its names and their order
        nt = barycentric((g=2.0,h=3.0,i=5.0))
        @test nt isa NamedTuple{(:g,:h,:i)}
        @test collect(nt) ≈ [0.2,0.3,0.5]
        @test collect(barycentric((i=5.0,g=2.0,h=3.0))) ≈ [0.5,0.2,0.3]
        @test sum(barycentric((a=1.0,b=1.0),4)) ≈ 4
        ## the tuple and NamedTuple methods agree
        @test collect(barycentric((g=2.0,h=3.0,i=5.0),7)) ≈ collect(barycentric((2.0,3.0,5.0),7))
        ## it is a projection: invariant to rescaling, and idempotent
        @test collect(barycentric(1000 .* x)) ≈ collect(b)
        @test collect(barycentric(b)) ≈ collect(b)
        ## a point already on the simplex is unchanged; one coordinate goes to n
        @test collect(barycentric((0.2,0.3,0.5))) ≈ [0.2,0.3,0.5]
        @test collect(barycentric((3.0,))) ≈ [1.0]
        @test collect(barycentric((a=3.0,),5)) ≈ [5.0]
        ## integer coordinates
        @test collect(barycentric((1,1,2))) ≈ [0.25,0.25,0.5]
    end

end
