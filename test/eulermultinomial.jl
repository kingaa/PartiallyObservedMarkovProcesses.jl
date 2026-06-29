using PartiallyObservedMarkovProcesses
using Statistics: mean, std
using Test

@info h1("testing Euler-multinomial facilities")

@testset verbose=true "eulermultinomial" begin

    d = EulerMultinomial(10,[1,2,3],0.1)
    @test d isa EulerMultinomial{3}

    d = EulerMultinomial(100,[1.0,2.0,3],0.2)
    @test d isa EulerMultinomial{3}
    x = rand(d)
    @test x isa Vector{<:Integer} && length(x)==3
    rand!(x,EulerMultinomial(100,[1.0,2.0,0],0.2))
    x = rand(EulerMultinomial(100,[0,0,0],0.1))
    @test all(x.==0)
    x = rand(EulerMultinomial(0,[1.0,2.0,3.0],0.1))
    @test all(x.==0)

    d = EulerMultinomial(10000,[1.0,2.0,3.0,0.0],0.3)
    @freeze 39586685 y = reduce(
        hcat,
        [rand(d) for _ ∈ 1:1000]
    )

    ma = mean(sum(y,dims=1))
    sa = std(sum(y,dims=1))
    @test ma-2*sa < 8347 < ma+2*sa
    m2 = mean(y[2,:]./y[1,:])
    s2 = std(y[2,:]./y[1,:])
    @test m2-2*s2 < 2 < m2+2*s2
    m3 = mean(y[3,:]./y[1,:])
    s3 = std(y[3,:]./y[1,:])
    @test m3-2*s3 < 3 < m3+2*s3
    @test sum(abs.(y[4,:]))==0

    @test_throws "invalid size" EulerMultinomial(-5,[1,2,3],0.1)
    @test_throws "invalid dt" EulerMultinomial(500,[1,2,3],-0.1)
    @test_throws "invalid rate[3]" EulerMultinomial(500,[1,2,-3],0.1)
    @test_throws "size mismatch" rand!(x,EulerMultinomial(50,[3.0,2.0],0.1))

end
