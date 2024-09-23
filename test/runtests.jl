using POMP
using Test

@testset "POMP.jl" begin
    @time include("test1.jl")
    @time include("test1.jl")
    @time include("gompertz.jl")
    @time include("gompertz.jl")
end
