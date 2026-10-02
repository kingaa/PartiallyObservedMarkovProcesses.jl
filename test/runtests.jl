using PartiallyObservedMarkovProcesses
using Test
using Crayons

h1 = crayon"bold blue"
h2 = s -> crayon"!bold light_yellow"("- "*s)

@testset verbose=true "POMP.jl" begin
    # The plotting tests must run before R draws any graphics (the
    # R-comparison tests below): R's graphics devices load the system
    # Glib, after which CairoMakie's own Glib fails to load.
    include("plots.jl")
    include("basic.jl")
    include("errors.jl")
    include("val_array.jl")
    include("helpers.jl")
    include("transformations.jl")
    include("eulermultinomial.jl")
    include("bake.jl")
    include("melt.jl")
    include("gompertz.jl")
    include("gompertz_kalman.jl")
    include("brown.jl")
    include("sir.jl")
    include("rmca.jl")
    include("drmca.jl")
    include("flow.jl")
    include("trajmatch.jl")
    include("pfilter.jl")
    include("iid.jl")
    include("weighted.jl")
    include("perturbn.jl")
    include("mif.jl")
    include("perturbn_support.jl")
    include("init_state.jl")
    include("core_fixes.jl")
    include("design.jl")
    include("slice.jl")
    include("mcap.jl")
    include("profile.jl")
    include("monitor.jl")
    include("speed1.jl")
end
