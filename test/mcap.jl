using PartiallyObservedMarkovProcesses
import PartiallyObservedMarkovProcesses as POMP
using DataFrames
using RCall
using Random
using Test

@info h1("mcap tests")

@testset verbose=true "mcap" begin

    Random.seed!(1517486)
    par = collect(range(1.0,3.0,length=40))
    ll = -5 .* (par .- 2.1).^2 .+ 0.3 .* randn(40)
    grid = collect(range(minimum(par),maximum(par),length=1000))

    @testset "loess against R (direct surface)" begin
        fit = POMP.loess(par,ll;span=0.75)
        sm = fit.(grid)
        @rput par ll grid
        R"""
        f_direct <- loess(ll ~ par, span=0.75, control=loess.control(surface="direct"))
        sm_direct <- predict(f_direct, newdata=grid)
        f_default <- loess(ll ~ par, span=0.75)
        sm_default <- predict(f_default, newdata=grid)
        """
        @rget sm_direct sm_default
        @test maximum(abs.(sm .- sm_direct)) < 1e-10
        # R's default interpolating surface differs only slightly
        @test maximum(abs.(sm .- sm_default)) < 0.05
        # span > 1 uses every point with an inflated bandwidth
        fit2 = POMP.loess(par,ll;span=1.5)
        R"""
        f2 <- loess(ll ~ par, span=1.5, control=loess.control(surface="direct"))
        sm2 <- predict(f2, newdata=grid)
        """
        @rget sm2
        @test maximum(abs.(fit2.(grid) .- sm2)) < 1e-10
        @test_throws r"same length" POMP.loess(par,ll[1:5])
        @test_throws r"finite" POMP.loess(par,vcat(ll[1:end-1],NaN))
        @test_throws r"span\*n" POMP.loess(par,ll;span=0.05)
    end

    @testset "mcap against R pomp::mcap (direct surface)" begin
        # pomp::mcap with its loess call switched to the direct surface,
        # which is what the Julia port computes; everything else verbatim.
        R"""
        library(pomp)
        mcap_direct <- function (logLik, parameter, level = 0.95, span = 0.75, Ngrid = 1000) {
            smooth_fit <- loess(logLik ~ parameter, span = span, control=loess.control(surface="direct"))
            parameter_grid <- seq(min(parameter), max(parameter), length.out = Ngrid)
            smoothed_logLik <- predict(smooth_fit, newdata = parameter_grid)
            smooth_arg_max <- parameter_grid[which.max(smoothed_logLik)]
            dist <- abs(parameter - smooth_arg_max)
            included <- dist < sort(dist)[trunc(span * length(dist))]
            maxdist <- max(dist[included])
            weights <- numeric(length(parameter))
            weights[included] <- (1 - (dist[included]/maxdist)^3)^3
            quadratic_fit <- lm(logLik ~ a + b, weights = weights,
                data = data.frame(logLik = logLik, b = parameter, a = -parameter^2))
            b <- unname(coef(quadratic_fit)["b"]); a <- unname(coef(quadratic_fit)["a"])
            m <- vcov(quadratic_fit)
            var_b <- m["b", "b"]; var_a <- m["a", "a"]; cov_ab <- m["a", "b"]
            se_mc_squared <- (1/(4 * a * a)) * (var_b - (2 * b/a) * cov_ab + (b * b/a/a) * var_a)
            se_stat_squared <- 1/2/a
            delta <- qchisq(level, df = 1) * (a * se_mc_squared + 0.5)
            logLik_diff <- max(smoothed_logLik) - smoothed_logLik
            ci <- range(parameter_grid[logLik_diff < delta])
            list(mle = smooth_arg_max, ci = ci, delta = delta,
                 se_stat = sqrt(se_stat_squared), se_mc = sqrt(se_mc_squared),
                 se = sqrt(se_mc_squared + se_stat_squared),
                 quadratic_max = b/(2*a), a = a, b = b,
                 c = unname(coef(quadratic_fit)[1]),
                 smoothed = smoothed_logLik)
        }
        md <- mcap_direct(ll, par)
        mp <- pomp::mcap(ll, par)
        """
        @rget md mp
        m = mcap(ll,par)
        @test m isa MCAP
        @test m.mle ≈ md[:mle] atol=1e-12
        @test m.ci[1] ≈ md[:ci][1] atol=1e-12
        @test m.ci[2] ≈ md[:ci][2] atol=1e-12
        @test m.delta ≈ md[:delta] rtol=1e-10
        @test m.se_stat ≈ md[:se_stat] rtol=1e-10
        @test m.se_mc ≈ md[:se_mc] rtol=1e-8
        @test m.se ≈ md[:se] rtol=1e-10
        @test m.quadratic_max ≈ md[:quadratic_max] rtol=1e-10
        @test m.coefs.a ≈ md[:a] rtol=1e-10
        @test m.coefs.b ≈ md[:b] rtol=1e-10
        @test m.coefs.c ≈ md[:c] rtol=1e-10
        @test maximum(abs.(m.fit.smoothed .- md[:smoothed])) < 1e-10
        @test m.fit.parameter == grid
        @test m.fit.quadratic ≈ m.coefs.c .+ m.coefs.b .* grid .- m.coefs.a .* grid.^2
        @test m.level == 0.95 && m.span == 0.75
        @test m.logLik == ll && m.parameter == par
        # and close to pomp::mcap as shipped (interpolating surface)
        step = grid[2]-grid[1]
        @test abs(m.mle - mp[:mle]) ≤ 10*step
        @test abs(m.ci[1] - mp[:ci][1]) ≤ 10*step
        @test abs(m.ci[2] - mp[:ci][2]) ≤ 10*step
        @test m.se_stat ≈ mp[:se_stat] rtol=0.01
        @test m.se_mc ≈ mp[:se_mc] rtol=0.01
        @test m.delta ≈ mp[:delta] rtol=0.01
        # the interval brackets the truth and the estimate
        @test m.ci[1] < 2.1 < m.ci[2]
        @test m.ci[1] ≤ m.mle ≤ m.ci[2]
        @test m.se ≈ sqrt(m.se_stat^2+m.se_mc^2)
        @test occursin("mle=",sprint(show,m))
    end

    @testset "mcap options and errors" begin
        m1 = mcap(ll,par;level=0.9,span=0.5,Ngrid=200)
        @test m1.level == 0.9 && nrow(m1.fit) == 200
        m2 = mcap(ll,par;level=0.99)
        @test m2.ci[2]-m2.ci[1] > mcap(ll,par;level=0.9).ci[2]-mcap(ll,par;level=0.9).ci[1]
        # integer inputs are accepted
        @test mcap(round.(Int,10 .* ll),round.(Int,10 .* par);span=1.0) isa MCAP
        @test_throws r"same length" mcap(ll,par[1:10])
        @test_throws r"finite" mcap(vcat(ll[1:end-1],-Inf),par)
        @test_throws r"level" mcap(ll,par;level=1.0)
        @test_throws r"Ngrid" mcap(ll,par;Ngrid=1)
        @test_throws r"span\*length" mcap(ll,par;span=0.01)
        # the quadratic fitting window needs span ≤ 1
        @test_throws r"\(0,1\]" mcap(ll,par;span=1.5)
        # too few points in the quadratic window: NaN results, with a warning
        few = @test_logs (:warn,r"carry weight") match_mode=:any mcap(ll[1:5],par[1:5];span=1.0)
        @test isnan(few.se) && all(isnan,few.ci)
        @test isfinite(few.mle)
        # a convex set of points gives no standard errors, with a warning
        mc = @test_logs (:warn,r"not concave") mcap(-ll,par)
        @test isnan(mc.se_stat) && isnan(mc.se)
        @test all(isnan,mc.ci)
        # a nonconcave fit whose total variance nonetheless comes out
        # positive (R reports se ≈ 0.171 and a negative cutoff here):
        # every uncertainty output is NaN, the fit itself is kept
        x = collect(range(1.0,3.0,length=40))
        mn = @test_logs (:warn,r"not concave") match_mode=:any mcap((x .- 2.1).^2 .+ 10 .* sin.(1:40),x)
        @test mn.coefs.a < 0 && isfinite(mn.mle)
        @test isnan(mn.se_stat) && isnan(mn.se_mc) && isnan(mn.se)
        @test isnan(mn.delta) && isnan(mn.quadratic_max) && all(isnan,mn.ci)
    end

    @testset "fits the data do not determine" begin
        # a profile over 5 values with 5 starts each: at the default span
        # only 2 distinct values carry weight in the quadratic window, so
        # the quadratic is not determined (R's lm gives NA)
        Random.seed!(11)
        x = repeat(collect(range(1.0,2.0,length=5)),inner=5)
        y = -6 .* (x .- 1.45).^2 .+ 0.3 .* randn(25)
        m = @test_logs (:warn,r"not determined") mcap(y,x)
        @test isfinite(m.mle)
        @test all(isnan,(m.coefs.c,m.coefs.a,m.coefs.b,m.quadratic_max,m.se_stat,m.se_mc,m.se,m.delta))
        @test all(isnan,m.ci) && all(isnan,m.fit.quadratic)
        # 8 points, span 0.5: a window of 3 points, 2 of them weighted
        x8 = collect(1.0:8.0)
        y8 = -(x8 .- 4.6).^2 .+ 0.3 .* sin.(1:8)
        m8 = @test_logs (:warn,r"not determined") mcap(y8,x8;span=0.5)
        @test isnan(m8.quadratic_max) && isnan(m8.se_stat) && isnan(m8.coefs.a)
        # the smoother's local quadratic needs 3 distinct values too: with
        # 4 values × 5 starts at the default span, the windows near the ends
        # hold 2 (R falls back on a pseudoinverse); this is refused
        x4 = repeat([1.0,1.5,2.0,2.5],inner=5)
        y4 = -4 .* (x4 .- 1.6).^2 .+ 0.5 .* randn(20)
        @test_throws r"`span` = 0.75 is too small" mcap(y4,x4)
        # with q = span*n points in a window, the farthest has weight zero:
        # q = 3 can never determine a local quadratic
        @test_throws r"`span\*n` must be at least 4" POMP.loess(collect(1.0:10.0),sin.(1:10);span=0.3)
        # 6 values × 3 starts at span 1: every window holds at least 3
        # distinct values, and the fit is determined
        x6 = repeat(collect(range(1.0,2.5,length=6)),inner=3)
        y6 = -4 .* (x6 .- 1.6).^2 .+ 0.5 .* randn(18)
        m6 = mcap(y6,x6;span=1.0)
        @test isfinite(m6.quadratic_max) && isfinite(m6.se) && all(isfinite,m6.ci)
    end

    @testset "mcap does not depend on the units of the parameter" begin
        x = collect(range(1.0,3.0,length=40))
        y = -5 .* (x .- 2.1).^2 .+ 0.3 .* sin.(1:40)
        m1 = mcap(y,x)
        for c ∈ (1e-8,1e8)
            mc = mcap(y,c .* x)
            @test collect(mc.ci) ./ c ≈ collect(m1.ci) rtol=1e-8
            @test mc.mle/c ≈ m1.mle rtol=1e-8
            @test mc.se/c ≈ m1.se rtol=1e-6
            @test mc.se_mc/c ≈ m1.se_mc rtol=1e-6
            @test mc.delta ≈ m1.delta rtol=1e-6
            @test mc.coefs.a*c^2 ≈ m1.coefs.a rtol=1e-6
        end
        # nor on its origin
        for d ∈ (1e6,-37.5)
            md = mcap(y,x .+ d)
            @test md.mle - d ≈ m1.mle atol=1e-8
            @test md.se ≈ m1.se rtol=1e-6
            @test md.delta ≈ m1.delta rtol=1e-6
            @test collect(md.ci) .- d ≈ collect(m1.ci) atol=1e-8
        end
        # nor on the order of the points
        p = randperm(MersenneTwister(1),40)
        mp = mcap(y[p],x[p])
        @test mp.mle == m1.mle && mp.ci == m1.ci && mp.fit.parameter == m1.fit.parameter
    end

end
