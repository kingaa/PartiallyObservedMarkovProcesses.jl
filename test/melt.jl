using POMP
using DataFrames
using Test

@testset "melt" begin
    d = melt(nothing);
    @test isa(d,DataFrame)
    @test size(d)==(0,0)
    d = melt((x=3,y=99));
    @test isa(d,DataFrame)
    @test size(d)==(1,2)
    @test names(d)==["x","y"]
    d = melt(fill((a=10,b=0,c=33),5,2));
    @test isa(d,DataFrame)
    @test size(d)==(10,5)
    @test names(d)==["id1","id2","a","b","c"]
    d = melt(fill((a=10,b=0,c=33),5,2),:i,:j);
    @test isa(d,DataFrame)
    @test size(d)==(10,5)
    @test names(d)==["i","j","a","b","c"]
end
