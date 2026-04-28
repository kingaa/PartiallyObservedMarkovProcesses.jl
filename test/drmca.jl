using PartiallyObservedMarkovProcesses
using PartiallyObservedMarkovProcesses.Examples
import PartiallyObservedMarkovProcesses as POMP
using RCall
using Test
using BenchmarkTools
using Random: seed!

@info h1("deterministic Rosenzweig-MacArthur model tests")

@testset "deterministic Rosenzweig-MacArthur model" begin

    seed!(875002133)

    @info h2("POMP.jl simulations (det Rosenzweig-MacArthur)")
    P = drmca()
    @test P isa POMP.PompObject
    @btime drmca()

    R"""
library(tidyverse,warn.conflicts=FALSE)
$(melt(P)) |>
   ggplot(aes(x=n,y=p))+
   geom_path()+
   scale_x_log10()+
   scale_y_log10()+
   theme_bw()
"""
    R"""ggsave(filename="drmca-01.png",width=7,height=7)"""

    R"""
$(melt(P)) |>
   select(time,n,p) |>
   pivot_longer(-time) |>
   ggplot(aes(x=time,y=value))+
   geom_line()+
   facet_grid(name~.,scales="free_y")+
   theme_bw()
"""
    R"""ggsave(filename="drmca-02.png",width=7,height=4)"""

    @info h2("POMP.jl simulation scaling (det Rosenzweig-MacArthur)")
    S = simulate(P,rmeasure=nothing,nsim=1)[1]
    @btime simulate_array($S,nsim=1)
    @btime simulate_array($S,nsim=10)
    @btime simulate_array($S,nsim=100)

end
