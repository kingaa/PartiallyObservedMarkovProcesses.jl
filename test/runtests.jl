using PartiallyObservedMarkovProcesses
using Test
using Crayons

h1 = crayon"bold blue"
h2 = s -> crayon"!bold light_yellow"("- "*s)

@testset verbose=true "POMP.jl" begin
    include("basic.jl")
    include("errors.jl")
    include("val_array.jl")
    include("helpers.jl")
    include("melt.jl")
    include("gompertz.jl")
    include("brown.jl")
    include("sir.jl")
    include("rmca.jl")
    include("drmca.jl")
    include("flow.jl")
    include("pfilter.jl")
    include("speed1.jl")
end
