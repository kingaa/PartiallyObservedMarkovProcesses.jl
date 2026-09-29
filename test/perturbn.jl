using PartiallyObservedMarkovProcesses
using Random: seed!
using DataFrames
using Chain: @chain
using Distributions: LogNormal
using Statistics: mean, std
using Test

@info h1("@perturbn tests")

@testset verbose=true "@perturbn" begin

    seed!(450237782)

    f = @perturbn a~LogNormal(0.1) b~ivp(LogNormal(1),0) p~LogitNormal(0.1) c~Normal(10) d~ivp(Normal(10),1) (e,f)~ivp(LogBaryNormal(1))

    x = f(1, 0, a=1, b=10, c=0, d=3, p=0.8, e=1, f=3, g=3, h=1, i=1)
    @test x.b != 10 && x.d == 3 && x.e+x.f ≈ 1
    x = f(1, 1, p=0.8, a=1, b=10, c=0, d=7, e=1, f=3, g=3, h=1, i=1)
    @test x.b == 10 && x.d != 7 && x.e == 1 && x.f == 3
    x = f(1, 10, a=1, b=10, c=0, p=0.8, d=5, e=1, f=3, g=3, h=1, i=1)
    @test x.b == 10 && x.d == 5

    @test_throws "Unrecognized perturbation specification" eval(:(@perturbn(p ~ LogCabin(0.1))))
    @test_throws "Unrecognized perturbation specification" eval(:(@perturbn(p ~ ivp(LogCabin(0.1)))))

    @test_throws "proper specification is" eval(:(@perturbn 3~LogNormal(1)))
    @test_throws "proper specification is" eval(:(@perturbn (a,b)~LogitNormal(1)))
    @test_throws "proper specification is" eval(:(@perturbn (a,)~LogBaryNormal(1)))
    @test_throws "proper `@perturbn` specification" eval(:(@perturbn x))
    @test_throws "proper `@perturbn` specification" eval(:(@perturbn x=LogNormal(5)))

    f = @perturbn(
        a~LogNormal(0.1),
        b~LogitNormal(0.1),
        (c,d,e)~LogBaryNormal(0.1),
        f~Normal(0),
        g~Normal(1),
        h ~ ivp(LogNormal(0.1),1),
        (i,j)~ivp(LogBaryNormal(0.1))
    )

    expand_grid(; kwargs...) = begin
        names, vals = keys(kwargs), values(kwargs)
        DataFrame(NamedTuple{names}(t) for t in Iterators.product(vals...))
    end

    nrep = 10000
    rtol = 0.03

    seed!(1138558856)

    tests = @chain expand_grid(
        scale=[1,2], lag=[0,1,5], rep=1:nrep,
        a=100, b=0.3, c=1, d=1, e=8, f=200, g=88, h=23, i=0.3, j=0.7
    ) begin
        transform(AsTable(Not([:scale,:lag])) => ByRow(identity) => :data)
        transform([:scale,:lag,:data] => ByRow((scale,lag,data,) -> f(scale,lag;data...)) => :data)
        transform(:data => AsTable)
        select(Not(:data,:rep))
        transform(
            AsTable([:c,:d,:e]) => ByRow(sum) => :Σcde,
            AsTable([:i,:j]) => ByRow(sum) => :Σij,
        )
        groupby([:lag,:scale])
        combine(
            :a => (x -> exp(mean(log.(x)))/100) => :μa,
            [:a,:scale] => ((x,s) -> std(log.(x))/0.1/mean(s)) => :σa,
            :b => (x -> expit(mean(logit.(x)))/0.3) => :μb,
            [:b,:scale] => ((x,s) -> std(logit.(x))/0.1/mean(s)) => :σb,
            :Σcde => mean => :μΣcde,
            :Σcde => std => :σΣcde,
            :d => (x -> exp(mean(log.(x)))/0.1) => :μd,
            :e => (x -> exp(mean(log.(x)))/0.8) => :μe,
            :f => (x -> mean(x)/200) => :μf,
            :f => std => :σf,
            :g => (x -> mean(x)/88) => :μg,
            [:g,:scale] => ((x,s) -> std(x)/mean(s)) => :σg,
            :h => (x -> exp(mean(log.(x)))/23) => :μh,
            [:h,:scale] => ((x,s) -> std(log.(x))/0.1/mean(s)) => :σh,
            :Σij => mean => :μΣij,
            :Σij => std => :σΣij,
            :i => (x -> exp(mean(log.(x)))/0.3) => :μi,
            [:i,:scale] => ((x,s) -> std(log.(x))/0.1/mean(s)) => :σi,
        )
        transform(
            :μa => ByRow(x -> isapprox(x, 1.0; rtol)) => :μa,
            :σa => ByRow(x -> isapprox(x, 1.0; rtol)) => :σa,
            :μb => ByRow(x -> isapprox(x, 1.0; rtol)) => :μb,
            :σb => ByRow(x -> isapprox(x, 1.0; rtol)) => :σb,
            :μΣcde => ByRow(x -> isapprox(x, 1.0; rtol)) => :μΣcde,
            :σΣcde => ByRow(x -> isapprox(x+1, 1.0; rtol)) => :σΣcde,
            :μd => ByRow(x -> isapprox(x, 1.0; rtol)) => :μd,
            :μe => ByRow(x -> isapprox(x, 1.0; rtol)) => :μe,
            :μf => ByRow(==(1.0)) => :μf,
            :σf => ByRow(==(0.0)) => :σf,
            :μg => ByRow(x -> isapprox(x, 1.0; rtol)) => :μg,
            :σg => ByRow(x -> isapprox(x, 1.0; rtol)) => :σg,
            :μh => ByRow(x -> isapprox(x, 1.0; rtol)) => :μh,
            :μi => ByRow(x -> isapprox(x, 1.0; rtol)) => :μi,
            [:σh, :lag] =>
                ByRow((s,lag) -> (lag != 1 || isapprox(s,1; rtol)) && (lag == 1 || isapprox(s+1,1; rtol))) =>
                :σh,
            [:σi, :lag] =>
                ByRow((s,lag) -> (lag != 0 || isapprox(s,1; rtol)) && (lag == 0 || isapprox(s+1,1; rtol))) =>
                :σi,
            :μΣij => ByRow(x -> isapprox(x, 1.0; rtol)) => :μΣij,
            :σΣij => ByRow(x -> isapprox(x+1, 1.0; rtol)) => :σΣij,
        )
        stack(Not(:lag,:scale),value_name=:pass)
    end

    @test all(tests.pass)
    
end
