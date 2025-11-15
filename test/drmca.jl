using PartiallyObservedMarkovProcesses
using PartiallyObservedMarkovProcesses.Examples
import PartiallyObservedMarkovProcesses as POMP
using RCall
using Test
using Random: seed!

println("- deterministic Rosenzweig-MacArthur model tests")

@testset "deterministic Rosenzweig-MacArthur model" begin

    seed!(875002133)

    println("    POMP.jl simulations (det Rosenzweig-MacArthur)")
    P = drmca()
    @time P = drmca()
    @time P = drmca()
    @time P = drmca()
    @test isa(P,POMP.PompObject)
    println("    ",P)

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

    println("    POMP.jl simulation scaling (det Rosenzweig-MacArthur)")
    S = simulate(P,rmeasure=nothing,nsim=1)[1]
    @time S1 = simulate(S,nsim=1)
    @time S1 = simulate(S,nsim=1)
    @time S1 = simulate(S,nsim=10)
    @time S1 = simulate(S,nsim=100)
    @time S1 = simulate(S,nsim=1000)

end
