using POMP
using Test

@testset "POMP.jl" begin
    include("test1.jl")
    include("gompertz.jl")
    include("pfilter.jl")
    include("speed1.jl")
end
