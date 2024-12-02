using POMP
using RCall
using Test
using Random: seed!

seed!(875002133)

@testset "Rosenzweig-MacArthur model" begin

    println("POMP.jl simulations (Rosenzweig-MacArthur)")
    @time P = rmca()
    @time P = rmca()
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

    R"""ggsave(filename="rmca-01.png",width=7,height=7)"""

    R"""
$(melt(P)) |>
   select(time,n,p) |>
   pivot_longer(-time) |>
   ggplot(aes(x=time,y=value))+
   geom_line()+
   facet_grid(name~.,scales="free_y")+
   theme_bw()
"""

    R"""ggsave(filename="rmca-02.png",width=7,height=4)"""

    println("POMP.jl pfilter (Rosenzweig-MacArthur)")
    P = rmca(δt=0.1,σ=0.1,times=range(1.0,20.0,step=1.0))
    @time Pf = pfilter(P,Np=1000)
    @time Pf = pfilter(P,Np=1000)
    @test isa(Pf,POMP.PfilterdPompObject)
    println(
        "POMP.jl likelihood estimate (Rosenzweig-MacArthur): ",
        round(Pf.logLik,digits=2)
    )

end
