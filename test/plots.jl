using PartiallyObservedMarkovProcesses
import PartiallyObservedMarkovProcesses as POMP
using PartiallyObservedMarkovProcesses.Examples: gompertz
using DataFrames
using Distributions
using AlgebraOfGraphics
using CairoMakie
import Random
using Test

@info h1("plots (AlgebraOfGraphics extension)")

@testset verbose=true "plots" begin

    Random.seed!(1953273041)

    rin = function(;x0,_...)
        (x=rand(Poisson(x0)),)
    end
    rlin = function (;t,a,x,_...)
        (x=rand(Poisson(a*x)),)
    end
    rmeas = function (;x,k,_...)
        (y=rand(NegativeBinomial(k,k/(k+x))),)
    end
    logdmeas = function (;x,y,k,_...)
        logpdf(NegativeBinomial(k,k/(k+x)),y)
    end

    p1 = (a=1.5,k=7.0,x0=5.0)

    P = simulate(
        t0=0,
        times=0:20,
        params=p1,
        rinit=rin,
        rprocess=discrete_time(rlin,dt=1),
        rmeasure=rmeas,
        logdmeasure=logdmeas
    )[1]

    ptb = @perturbn(a ~ LogNormal(0.02), k ~ LogNormal(0.05), x0 ~ ivp(LogNormal(0.1)))
    mf1 = mif(P;Nmif=5,Np=50,perturbations=ptb,cooling=geometric_cooling(0.5))
    mf2 = mif(P;Nmif=5,Np=50,perturbations=ptb,cooling=hyperbolic_cooling(0.5))
    mon2 = monitor(mf2;Np=50,seed=1,every=2)

    ext = Base.get_extension(PartiallyObservedMarkovProcesses,:PartiallyObservedMarkovProcessesAoGExt)
    @test !isnothing(ext)

    # what a panel of a drawn figure actually shows: its plots of a given
    # kind (:scatter, :lines, :errorbars, :vlines, :hlines), by their data
    shown(fg, i, kind) = [p[1][] for p ∈ fg.grid[i].axis.scene.plots if CairoMakie.Makie.plotkey(p) == kind]
    pts(x, y) = [CairoMakie.Point2(a,b) for (a,b) ∈ zip(x,y)]
    plotsof(fg, i, kind) = [p for p ∈ fg.grid[i].axis.scene.plots if CairoMakie.Makie.plotkey(p) == kind]

    @testset "traceplot" begin
        d, order = ext.trace_data([mf1],nothing,nothing)
        @test order == ["logLik","a","k","x0"]
        @test Set(d.variable) == Set(order)
        # the last trace row has no mif log likelihood
        @test count(==("logLik"),d.variable) == 5
        d, order = ext.trace_data([mf2],(:k,),[mon2])
        @test order == ["logLik","monitor logLik","k"]
        # the monitor shares the iteration axis of the traces
        @test d.iteration[d.variable .== "monitor logLik"] == [1,3,5,6]
        @test d.iteration[d.variable .== "k"] == 1:6
        @test_throws r"not a parameter" ext.trace_data([mf1],(:bogus,),nothing)
        @test_throws r"no `mif` results" ext.trace_data(POMP.MifdPompObject[],nothing,nothing)
        @test_throws r"one `monitor`" ext.trace_data([mf1,mf2],nothing,[mon2])
        # a monitor table computed for other `mif` results is refused,
        # whether they are a continuation or another run of equal length
        mf1c = mif(mf1;Nmif=2)
        @test_throws r"one `monitor`" ext.trace_data([mf1,mf1c],nothing,[monitor([mf1,mf1c];Np=20,seed=1)])
        @test_throws r"computed for other" ext.trace_data([mf1c],nothing,[monitor([mf1,mf1c];Np=20,seed=1)])
        @test_throws r"computed for other" ext.trace_data([mf1],nothing,[monitor(mf2;Np=20,seed=1)])
        # by default, the parameters perturbed in any of the runs
        ma = mif(P;Nmif=2,Np=20,perturbations=@perturbn(a ~ LogNormal(0.02)),cooling=geometric_cooling(0.5))
        mk = mif(P;Nmif=2,Np=20,perturbations=@perturbn(k ~ LogNormal(0.05)),cooling=geometric_cooling(0.5))
        @test ext.trace_data([ma,mk],nothing,nothing)[2] == ["logLik","a","k"]
        fg = traceplot(mf1)
        @test fg isa AlgebraOfGraphics.FigureGrid
        @test size(fg.grid) == (2,2)
        # the first panel is the mif log likelihood against the iteration
        t = traces(mf1)
        @test shown(fg,1,:lines) == [pts(t.iteration[1:end-1],t.logLik[1:end-1])]
        save("mif-01.png",fg)
        @test isfile("mif-01.png")
        # a parameter named `iteration` would replace the iteration number,
        # the horizontal axis of every panel
        Mi = mif(P;params=merge(p1,(iteration=2.0,)),Nmif=1,Np=10,
            perturbations=@perturbn(a ~ LogNormal(0.02)),cooling=geometric_cooling(0.5))
        @test_throws r"`iteration` is reserved" traceplot(Mi)
        # a single parameter may be named on its own
        @test ext.trace_data([mf1],:k,nothing)[2] == ["logLik","k"]
        fg = traceplot([mf1,mf2];pars=(:a,:k))
        @test size(fg.grid) == (2,2)
        # each parameter's panel shows that parameter's trace, for each run
        t1, t2 = traces(mf1), traces(mf2)
        @test Set(shown(fg,CartesianIndex(1,2),:lines)) == Set([pts(t1.iteration,t1.a),pts(t2.iteration,t2.a)])
        @test Set(shown(fg,CartesianIndex(2,1),:lines)) == Set([pts(t1.iteration,t1.k),pts(t2.iteration,t2.k)])
        fg = traceplot(mf2;monitor=mon2)
        # the second panel is the monitor, at its own iterations
        @test shown(fg,CartesianIndex(1,2),:lines) == [pts(mon2.iteration,mon2.loglik)]
        save("mif-02.png",fg)
        @test isfile("mif-02.png")
        @test_throws r"result of `mif`" traceplot(1)
    end

    @testset "filterplot" begin
        # the data: one row per run, time, and panel
        d = ext.filter_data([mf1])
        ess = d[d.variable .== "effective sample size",:]
        cll = d[d.variable .== "conditional log likelihood",:]
        @test ess.time == Float64.(times(mf1)) && ess.value == eff_sample_size(mf1)
        @test cll.time == Float64.(times(mf1)) && cll.value == cond_logLik(mf1)
        # a conditional log likelihood of -Inf becomes a gap (NaN)
        pdeg = pfilter(pomp(P;logdmeasure=(;_...) -> -Inf);Np=10)
        dd = ext.filter_data([pdeg])
        @test all(isnan,dd.value[dd.variable .== "conditional log likelihood"])
        # what is drawn: for each run, the effective sample size (top) and
        # the conditional log likelihood (bottom) against time
        fg = filterplot([mf1,mf2])
        @test fg isa AlgebraOfGraphics.FigureGrid
        @test size(fg.grid) == (2,1)
        t = Float64.(times(mf1))
        @test Set(shown(fg,1,:lines)) == Set([pts(t,eff_sample_size(m)) for m ∈ (mf1,mf2)])
        @test Set(shown(fg,2,:lines)) == Set([pts(t,cond_logLik(m)) for m ∈ (mf1,mf2)])
        save("filter-01.png",fg)
        @test isfile("filter-01.png")
        # a `pfilter` result, alone
        pf = pfilter(P;Np=50)
        fp = filterplot(pf)
        @test shown(fp,1,:lines) == [pts(Float64.(times(pf)),eff_sample_size(pf))]
        @test shown(fp,2,:lines) == [pts(Float64.(times(pf)),cond_logLik(pf))]
        # the time axis is labelled with the model's time variable
        pg = pfilter(gompertz();Np=20,params=(r=4.5,K=210.0,σₚ=0.7,σₘ=0.1,X0=150.0))
        labels = [string(c.text[]) for c ∈ filterplot(pg).figure.content if c isa CairoMakie.Makie.Label]
        @test "year" ∈ labels
        @test_throws r"result of `pfilter` or `mif`" filterplot(1)
        @test_throws r"result of `pfilter` or `mif`" filterplot([1,2])
        @test_throws r"no results" filterplot(POMP.PfilterdPompObject[])
    end

    @testset "sliceplot" begin
        d = slice_design(p1; a=[1.2,1.5,1.8], k=[4.0,7.0,10.0])
        s = slice(P,d;Np=50,nreps=2)
        fg = sliceplot(s)
        @test fg isa AlgebraOfGraphics.FigureGrid
        @test size(fg.grid) == (1,2)
        save("slice-01.png",fg)
        @test isfile("slice-01.png")
        # the points and error bars drawn, on a small known slice: a
        # missing standard error gives no error bar
        s0 = DataFrame(a=[1.0,2.0,1.5,1.5],k=[7.0,7.0,5.0,9.0],slice=[:a,:a,:k,:k],
            loglik=[-10.0,-12.0,-11.0,-13.0],se=[0.5,0.25,NaN,1.0])
        f0 = sliceplot(s0)
        @test shown(f0,1,:scatter) == [pts([1.0,2.0],[-10.0,-12.0])]
        @test shown(f0,2,:scatter) == [pts([5.0,9.0],[-11.0,-13.0])]
        @test [[v[3] for v ∈ e] for e ∈ shown(f0,1,:errorbars)] == [[0.5,0.25]]
        @test [[(v[1],v[3]) for v ∈ e] for e ∈ shown(f0,2,:errorbars)] == [[(9.0,1.0)]]
        # without a standard error column there are no error bars
        f1 = sliceplot(s0[:,Not(:se)])
        @test size(f1.grid) == (1,2)
        @test isempty(shown(f1,1,:errorbars)) && isempty(shown(f1,2,:errorbars))
        # options are passed to `draw`
        @test sliceplot(s;axis=(width=200,height=150)) isa AlgebraOfGraphics.FigureGrid
        @test_throws r"`slice` and `loglik` columns" sliceplot(DataFrame(a=[1.0],loglik=[1.0]))
        # a row whose estimate is -Inf (a zero likelihood) does not stop the plot
        @test sliceplot(DataFrame(a=[1.0,2.0,1.5],slice=fill(:a,3),loglik=[-3.0,-Inf,-4.0])) isa AlgebraOfGraphics.FigureGrid
        @test_throws r"data frame" sliceplot(1)
    end

    @testset "mcapplot" begin
        par = collect(range(1.0,3.0,length=30))
        ll = -5 .* (par .- 2.1).^2 .+ 0.3 .* randn(30)
        m = mcap(ll,par)
        fg = mcapplot(m)
        @test fg isa AlgebraOfGraphics.FigureGrid
        # the profile points, the smoothed and quadratic fits, the
        # estimate, both ends of the interval, and the cutoff
        @test shown(fg,1,:scatter) == [pts(m.parameter,m.logLik)]
        fits = shown(fg,1,:lines)
        @test pts(m.fit.parameter,m.fit.smoothed) ∈ fits
        @test pts(m.fit.parameter,m.fit.quadratic) ∈ fits
        # the dashed line is the quadratic fit, the solid one the smoothed profile
        dashed = [p[1][] for p ∈ plotsof(fg,1,:lines) if !isnothing(p.linestyle[])]
        solid = [p[1][] for p ∈ plotsof(fg,1,:lines) if isnothing(p.linestyle[])]
        @test dashed == [pts(m.fit.parameter,m.fit.quadratic)]
        @test solid == [pts(m.fit.parameter,m.fit.smoothed)]
        @test sort(reduce(vcat,shown(fg,1,:vlines))) ≈ sort([m.mle,m.ci...])
        @test only(only(shown(fg,1,:hlines))) ≈ maximum(m.fit.smoothed)-m.delta
        # a profile without an interval draws neither interval nor cutoff
        mc = @test_logs (:warn,r"not concave") match_mode=:any mcap(-ll,par)
        fc = mcapplot(mc)
        @test reduce(vcat,shown(fc,1,:vlines)) == [mc.mle]
        @test isempty(shown(fc,1,:hlines))
        # an undetermined quadratic (all NaN) still draws the points, the
        # smoothed profile, and the estimate, with no interval
        xu = repeat(collect(range(1.0,2.0,length=5)),inner=5)
        yu = -6 .* (xu .- 1.45).^2 .+ 0.3 .* randn(Random.MersenneTwister(11),25)
        mu = @test_logs (:warn,r"not determined") mcap(yu,xu)
        fu = mcapplot(mu)
        @test all(isnan,mu.fit.quadratic)
        @test reduce(vcat,shown(fu,1,:vlines)) == [mu.mle] && isempty(shown(fu,1,:hlines))
        save("mcap-01.png",fg)
        @test isfile("mcap-01.png")
        @test_throws r"MCAP" mcapplot(ll)
    end

end
