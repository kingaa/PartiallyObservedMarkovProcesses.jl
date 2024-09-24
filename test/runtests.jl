using POMP
using Test

@testset "POMP.jl" begin
    include("test1.jl")
    include("gompertz.jl")
end
