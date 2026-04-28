using PartiallyObservedMarkovProcesses
using DataFrames
using Test

@info h1("testing `melt` methods")

@testset verbose=true "melt methods" begin
    d = melt(nothing);
    @test d isa DataFrame
    @test size(d)==(0,0)
    d = melt((x=3,y=99));
    @test d isa DataFrame
    @test size(d)==(1,2)
    @test names(d)==["x","y"]
    d = melt(fill((a=10,b=0,c=33),5,2));
    @test d isa DataFrame
    @test size(d)==(10,5)
    @test names(d)==["id1","id2","a","b","c"]
    d = melt(fill((a=10,b=0,c=33),5,2),:i,:j);
    @test d isa DataFrame
    @test size(d)==(10,5)
    @test names(d)==["i","j","a","b","c"]
    d = melt(rand(5,3,2),:i,:j,:k);
    @test d isa DataFrame
    @test size(d)==(30,4)
    @test names(d)==["i","j","k","value"]
end
