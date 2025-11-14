using PartiallyObservedMarkovProcesses
using PartiallyObservedMarkovProcesses.Examples
import PartiallyObservedMarkovProcesses as POMP
using RCall
using Test
using Random: seed!

@testset "deterministic Rosenzweig-MacArthur model" begin

    seed!(875002133)

    println("PartiallyObservedMarkovProcesses.jl simulations (det Rosenzweig-MacArthur)")
    P = drmca()
    @time P = drmca()
    @time P = drmca()
    @time P = drmca()
    @test isa(P,POMP.PompObject)

    R"""
library(tidyverse)
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

end
