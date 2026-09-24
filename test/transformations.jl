using PartiallyObservedMarkovProcesses
import PartiallyObservedMarkovProcesses as POMP
using Test

@info h1("parameter transformations")

@testset verbose=true "parameter transformations" begin

    p = 0.2
    @test expit(logit(p)) ≈ p
    x = 4.0
    @test logit(expit(x)) ≈ x

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

    x = (2.0,3.0,5.0)
    b = barycentric(Int64,x,20)
    @test b isa NTuple{3,Int64}
    @test collect(b) ≈ [4,6,10]
    @test sum(b) == 20
    ## onto the simplex of total n
    @test collect(barycentric(Float64,x,10)) ≈ [2.0,3.0,5.0]
    @test collect(barycentric(Int64,(0.9,0.01,0.1),10000)) == [8911,99,990]
    @test sum(barycentric(Float64,(0.9,0.01,0.1),10000)) ≈ 10000
    ## a NamedTuple keeps its names and their order
    nt = barycentric(Int64,(g=2.0,h=3.0,i=5.0),10)
    @test nt isa NamedTuple{(:g,:h,:i)}
    @test collect(nt) ≈ [2,3,5]
    @test collect(barycentric(Int64,(i=5.0,g=2.0,h=3.0),10)) ≈ [5,2,3]
    @test sum(barycentric(Int64,(a=1.0,b=1.0),4)) == 4
    ## the tuple and NamedTuple methods agree
    @test collect(barycentric(Int64,(g=2.0,h=3.0,i=5.0),7)) == collect(barycentric(Int64,(2.0,3.0,5.0),7))
    ## it is a projection: invariant to rescaling, and idempotent
    @test collect(barycentric(Int64, 1000 .* x, 20)) == collect(b)
    @test collect(barycentric(Int64,b,20)) ≈ collect(b)

end
