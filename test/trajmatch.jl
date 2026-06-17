using PartiallyObservedMarkovProcesses
using PartiallyObservedMarkovProcesses.Examples
import PartiallyObservedMarkovProcesses as POMP
import OrdinaryDiffEqLowOrderRK: Euler
using Test
using Random: seed!

@info h1("trajectory matching tests")

@testset "trajectory matching" begin

    seed!(875002133)

    @test_throws "Incorrect call" traj_match_objfun("test")

    P = drmca()
    @test P isa POMP.PompObject

    f1 = traj_match_objfun(
        P,
        rmeasure=function(;X,Y,_...)
            (n=exp(X),p=exp(Y),)
        end
    )
    @test f1 isa Function
    @test round(f1(coef(P)),sigdigits=5)==4620.3
    @test_throws MethodError f1()
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
    @test_throws UndefKeywordError f3([1.0, 10000.0])

    f4 = traj_match_objfun(
        P,(:r,:σ),
        rmeasure=function(;X,Y,_...)
            (n=0.0,p=exp(Y),)
        end
    )
    @test f4([1.0, 2.0])==Inf

    f5 = traj_match_objfun(P,(:N0,:P0,:K,:A,:σ))
    @test_throws DomainError f5([-100.0, 2.0, 100, 0.0, 1.0])
    @test_throws MethodError f5(["a", 2.0, 100, 0.0, 1.0])
    @test @test_logs (:warn,r"unsuccessful vectorfield") match_mode=:any f5([ 100.0, 2.0, 100, 0.0, -1.0])==Inf

    f6 = traj_match_objfun(
        P, (:r,:K),
        rprocess=vectorfield(P.rprocess.vf,Euler(),dt=1)
    )
    @test @test_logs (:warn,r"instability") (:warn,r"unsuccessful vectorfield") f6([100.0,0.1])==Inf

    f7 = traj_match_objfun(
        P,(:r,:σ),
        rmeasure=function(;X,Y,_...)
            error("this is a test, this is only a test")
            (n=0.0,p=exp(Y),)
        end,
        whitelist=Union{ErrorException}
    )
    @test @test_logs (:error,r"whitelisted error") f7([100.0,0.1])==Inf

end
