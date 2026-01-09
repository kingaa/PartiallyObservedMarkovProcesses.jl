using PartiallyObservedMarkovProcesses
import PartiallyObservedMarkovProcesses as POMP
using Test

@info "testing val_array"

@testset "val_array" begin
    y = fill((a=7,b=99,c="bob"),25,1)
    @test POMP.val_array("yes")==["yes"]
    @test_throws "size mismatch" POMP.val_array(y,11,2)
    @test size(POMP.val_array(3,1,1,1))==(1,1,1,1)
    @test size(POMP.val_array(rand(3,2,2)))==(12,)
end
