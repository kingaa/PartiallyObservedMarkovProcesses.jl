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
    @test :x ∉ paramsymbs(P)
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

    @test_throws "Incorrect call" traj_match_objfun("test")

    f1 = traj_match_objfun(
        P,
        rmeasure=function(;X,Y,_...)
            (n=exp(X),p=exp(Y),)
        end
    )
    @test f1 isa Function
    @test round(f1(coef(P)),sigdigits=5)==4620.3
    @test_throws "no method matching" f1()
    @test_throws "incorrect argument length" f1([coef(P,:r,:K)...])

    f2 = traj_match_objfun(
        P,(:r,:σ),
        rmeasure=function(;X,Y,_...)
            (n=exp(X),p=exp(Y),)
        end
    )
    @test f2 isa Function
    @test f2((r=1.0,σ=2.0)) == f2([1.0, 2.0])
    @test f2([1.0, 2]) == f2([1, 2])
    @test f2(coef(P,:r,:σ))==f2(coef(P,:σ,:r))
    @test f2(coef(P,:r,:σ))!=f2([coef(P,:σ,:r)...])

    f3 = traj_match_objfun(
        P,(:r,:sigma),
        params=nothing,
        rmeasure=function(;X,Y,_...)
            (n=exp(X),p=exp(Y),)
        end
    )
    @test_throws r"keyword argument .* not assigned" f3([1.0, 10000.0])

    f4 = traj_match_objfun(
        P,(:r,:sigma),
        rmeasure=function(;X,Y,_...)
            @warn("rmeasure warning!")
            error("rmeasure error!")
            (n=exp(X),p=exp(Y),)
        end
    )
    val = @test_logs (:warn,"rmeasure warning!") (:warn,r"rmeasure error!") f4((1.0,10000.0))
    @test val==Inf

    f5 = traj_match_objfun(
        P,(:r,:sigma),
        rmeasure=function(;X,Y,_...)
            @info "rmeasure info"
            @warn "rmeasure warning!"
            (n=exp(X),p=exp(Y),)
        end
    )
    val = @test_logs (:warn,"rmeasure warning!") (:info,"rmeasure info") match_mode=:any f5((1.5,10000.0))
    @test isfinite(val)

    f6 = traj_match_objfun(
        P,(:r,:sigma),
        rmeasure=function(;X,Y,_...)
            (n=exp(X),p=1+"3",)
        end
    )
    val = @test_logs (:warn,r"MethodError") f6((1.5,10000.0))
    @test val==Inf

end
