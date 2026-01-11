using PartiallyObservedMarkovProcesses
using PartiallyObservedMarkovProcesses.Examples
using Test

@info "testing error capture"

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
        logdmeasure=logdmeas,
        logdprior=logdpri
    )

    x0 = rinit(P);
    x = rprocess(P,x0=x0);
    y = rmeasure(P,x=x);

    @test_throws "Incorrect call" simulate(P,rprocess=onestep("bob"))
    @test_throws "Incorrect call" simulate(P,rprocess=euler("bob"))
    @test_throws "Incorrect call" simulate(P,rprocess=vectorfield("bob"))
    @test_throws "Incorrect call" discrete_time("bob")
    @test_throws "Incorrect call" pomp("bob")
    @test_throws "Incorrect call" pfilter("bob")
    @test_throws "Incorrect call" simulate("bob")

    P1 = pomp(
        sir(),
        rinit=(x)->x,
        rprocess=onestep((x)->x),
        rmeasure=(x)->x,
        logdmeasure=(x)->1,
        logdprior=(x)->1
    );
    x0 = rinit(sir());
    x = rprocess(sir());
    @test_throws "no method matching" rinit(P1)
    @test_throws "no method matching" rinit!(P1,x0)
    @test_throws "no method matching" rprocess(P1,x0=x0)
    @test_throws "no method matching" rmeasure(P1,x=x)
    @test_throws "no method matching" logdmeasure(P1,x=x)
    @test_throws "no method matching" logdprior(P1)

end
