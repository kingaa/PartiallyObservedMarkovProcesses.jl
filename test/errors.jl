using POMP
using Test

@testset verbose=true "error messages" begin

    include("test1.jl")

    @test_throws "must be no later than" pomp(parus_data,times=:year,t0=1999)
    @test_throws "times must be nondecreasing" pomp(sort(parus_data,:pop),times=:year,t0=1940)
    @test_throws "must be of the same length" pomp([(a=1,);(a=2,);(c=3,)];times=[1;2],t0=0)

    P = pomp(
        times=1:21,
        t0=0,
        params=theta,
        rinit=rin,
        rprocess=discrete_time(rlin,dt=1),
        rmeasure=rmeas,
        logdmeasure=logdmeas
    )

    x0 = rinit(P);
    x = rprocess(P,x0=x0);
    y = rmeasure(P,x=x);

    P = pomp(P,rinit=function(;_...) error("yikes!") end);
    @test_throws "in `rinit`: yikes!" rinit(P,params=(x₀=3,))
    @test_throws "in `rinit!`: yikes!" rinit!(P,x0,params=(x₀=3,))
    P = pomp(P,rprocess=onestep(function(;_...) error("yikes!") end));
    @test_throws "in `rprocess!`: yikes!" rprocess(P,params=(x₀=3,),x0=x0)
    P = pomp(P,rmeasure=function(;_...) error("yikes!") end);
    @test_throws "in `rmeasure`: yikes!" rmeasure(P,params=(a=1,),x=x)
    P = pomp(P,logdmeasure=function(;_...) error("yikes!") end);
    @test_throws "in `logdmeasure!`: yikes!" logdmeasure(P,params=(a=1,),x=x,y=y)

    @test_throws "Incorrect call" simulate(P,rprocess=onestep("bob"))
    @test_throws "Incorrect call" simulate(P,rprocess=euler("bob"))
    @test_throws "Incorrect call" discrete_time("bob")
    @test_throws "Incorrect call" pomp("bob")    
    @test_throws "Incorrect call" pfilter("bob")
    @test_throws "Incorrect call" simulate("bob")    

end
