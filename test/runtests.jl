using POMP
using Test

@testset verbose=true "POMP.jl" begin
    include("test1.jl")
    include("gompertz.jl")
    include("sir.jl")
    include("rmca.jl")
    include("pfilter.jl")
    include("speed1.jl")
end
