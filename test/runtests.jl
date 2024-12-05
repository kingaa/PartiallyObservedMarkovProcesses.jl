using POMP
using Test

@testset verbose=true "POMP.jl" begin
    include("basic.jl")
    include("errors.jl")
    include("val_array.jl")
    include("helpers.jl")
    include("melt.jl")
    include("gompertz.jl")
    include("sir.jl")
    include("rmca.jl")
    include("pfilter.jl")
    include("speed1.jl")
end
