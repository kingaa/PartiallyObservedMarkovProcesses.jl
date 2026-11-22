using PartiallyObservedMarkovProcesses
import PartiallyObservedMarkovProcesses as POMP
using DataFrames
using Random
using Test

@info h1("design tests")

@testset verbose=true "designs" begin

    center = (a=1.5,k=7.0,x0=5.0)

    @testset "slice_design" begin
        d = slice_design(center; a=1.0:0.5:2.0, k=[3.0,9.0])
        @test d isa DataFrame
        @test size(d) == (5,4)
        @test propertynames(d) == [:a,:k,:x0,:slice]
        @test d.slice == [:a,:a,:a,:k,:k]
        @test d.a[1:3] == [1.0,1.5,2.0]
        @test all(d.k[1:3] .== 7.0)
        @test all(d.x0 .== 5.0)
        @test d.k[4:5] == [3.0,9.0]
        @test all(d.a[4:5] .== 1.5)
        @test eltype(d.a) == Float64
        # a scalar slice value is accepted
        @test nrow(slice_design(center; a=2.0)) == 1
        @test_throws r"does not appear in `center`" slice_design(center; bogus=1:3)
        @test_throws r"at least one slice" slice_design(center)
        @test_throws r"real numbers" slice_design((a=1.0,b="x"); a=1:2)
        # names used for output columns cannot be parameters
        @test_throws r"`se` is reserved" slice_design((a=1.0,se=2.0); a=1:2)
        @test_throws r"`slice` is reserved" slice_design((a=1.0,slice=2.0); a=1:2)
        # integer values stay integers; other values become Float64
        di = slice_design((N=50,p=0.3); p=[0.2,0.4])
        @test eltype(di.N) == Int && di.N == [50,50] && eltype(di.p) == Float64
        @test eltype(slice_design((N=50,p=0.3); N=[40,60]).N) == Int
        @test eltype(slice_design((N=50,p=0.3); N=[40.0,60.0]).N) == Int
        @test eltype(slice_design((N=50,p=0.3); N=[40.5]).N) == Float64
        @test eltype(slice_design((N=50,p=0.3); p=1).p) == Float64
        # values that the centre's integer type cannot hold exactly become
        # Float64, unchanged
        for (c,v) ∈ ((5,[1e20]),(5,[2.0^63]),(UInt8(5),[-1,5]),(UInt8(5),[300]),(true,[0,2]))
            dv = slice_design((N=c,p=0.3); N=v)
            @test eltype(dv.N) == Float64 && dv.N == Float64.(v)
        end
    end

    @testset "runif_design" begin
        lower = (k=3.0,x0=2.0)
        upper = (k=10.0,x0=8.0)
        d = runif_design(lower,upper,50;rng=MersenneTwister(1))
        @test size(d) == (50,2)
        @test propertynames(d) == [:k,:x0]
        @test all(3.0 .≤ d.k .≤ 10.0)
        @test all(2.0 .≤ d.x0 .≤ 8.0)
        d2 = runif_design(lower,upper,50;rng=MersenneTwister(1))
        @test d == d2
        # the draws themselves, column by column from the same stream
        rng = MersenneTwister(1)
        @test d.k == 3.0 .+ 7.0 .* rand(rng,50)
        @test d.x0 == 2.0 .+ 6.0 .* rand(rng,50)
        # names may be given in a different order in `upper`
        d3 = runif_design(lower,(x0=8.0,k=10.0),5;rng=MersenneTwister(2))
        @test propertynames(d3) == [:k,:x0]
        @test nrow(runif_design(lower,upper,0)) == 0
        @test_throws r"must match" runif_design(lower,(k=10.0,z=1.0),5)
        @test_throws r"at least as large" runif_design(lower,(k=1.0,x0=8.0),5)
        @test_throws r"nseq" runif_design(lower,upper,-1)
    end

    @testset "sobol_design" begin
        lower = (k=3.0,x0=2.0)
        upper = (k=10.0,x0=8.0)
        d = sobol_design(lower,upper,64)
        @test size(d) == (64,2)
        @test all(3.0 .≤ d.k .≤ 10.0)
        @test all(2.0 .≤ d.x0 .≤ 8.0)
        @test d == sobol_design(lower,upper,64)
        # `upper` is matched by name, not by position
        @test sobol_design(lower,(x0=8.0,k=10.0),64) == d
        # low discrepancy: every quarter of each range gets about a quarter
        # of the points (the sequence skips its initial zero point, so the
        # first 64 points are not an exact balanced block)
        for c ∈ (:k,:x0)
            lo, hi = lower[c], upper[c]
            q = [count(v -> lo+(j-1)*(hi-lo)/4 ≤ v < lo+j*(hi-lo)/4,d[!,c]) for j ∈ 1:4]
            @test all(14 .≤ q .≤ 18)
        end
        # jointly, not only margin by margin: each cell of a 4×4 grid on
        # the box holds 3 to 5 of the 64 points
        cell(v,lo,hi) = min(4,1+floor(Int,4*(v-lo)/(hi-lo)))
        C = zeros(Int,4,4)
        for r ∈ eachrow(d)
            C[cell(r.k,3.0,10.0),cell(r.x0,2.0,8.0)] += 1
        end
        @test all(3 .≤ C .≤ 5)
        # the points are those of the Sobol' sequence, scaled to the box
        s = POMP.SobolSeq(2)
        pts = [copy(POMP.next!(s,zeros(2))) for _ ∈ 1:64]
        @test d.k ≈ 3.0 .+ 7.0 .* first.(pts)
        @test d.x0 ≈ 2.0 .+ 6.0 .* last.(pts)
        @test nrow(sobol_design((a=0.0,),(a=1.0,),0)) == 0
    end

    @testset "profile_design" begin
        lower = (k=3.0,x0=2.0)
        upper = (k=10.0,x0=8.0)
        d = profile_design(a=[1.0,1.5,2.0];lower,upper,nprof=4,rng=MersenneTwister(3))
        @test size(d) == (12,3)
        @test propertynames(d) == [:a,:k,:x0]
        @test d.a == repeat([1.0,1.5,2.0],inner=4)
        @test all(3.0 .≤ d.k .≤ 10.0)
        @test metadata(d,"profiled") == [:a]
        # the starting values are the draws of `runif_design` from the same stream
        @test d.k == runif_design(lower,upper,12;rng=MersenneTwister(3)).k
        # two profiled variables: full grid, the first varying fastest
        d2 = profile_design(a=[1.0,2.0],b=[10.0,20.0,30.0];lower,upper,nprof=2,type=:sobol)
        @test size(d2) == (12,4)
        @test d2.a == repeat([1.0,2.0,1.0,2.0,1.0,2.0],inner=2)
        @test d2.b == repeat([10.0,10.0,20.0,20.0,30.0,30.0],inner=2)
        @test metadata(d2,"profiled") == [:a,:b]
        @test d2 == profile_design(a=[1.0,2.0],b=[10.0,20.0,30.0];lower,upper,nprof=2,type=:sobol)
        @test_throws r"at least one variable" profile_design(;lower,upper,nprof=2)
        @test_throws r"nprof" profile_design(a=[1.0];lower,upper,nprof=0)
        @test_throws r"type" profile_design(a=[1.0];lower,upper,nprof=1,type=:bogus)
        @test_throws r"must not appear in `lower`" profile_design(k=[1.0];lower,upper,nprof=1)
        @test_throws r"`slice` is reserved" profile_design(slice=[2.0,3.0];lower,upper,nprof=1)
        @test_throws r"`loglik` is reserved" profile_design(a=[1.0];lower=(loglik=0.0,),upper=(loglik=1.0,),nprof=1)
        # a grid written with integers is a grid of real values
        dn = profile_design(N=[40,60];lower,upper,nprof=2,rng=MersenneTwister(4))
        @test eltype(dn.N) == Float64 && dn.N == [40,40,60,60] && eltype(dn.k) == Float64
    end

end
