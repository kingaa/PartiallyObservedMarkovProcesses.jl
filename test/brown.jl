using PartiallyObservedMarkovProcesses
using PartiallyObservedMarkovProcesses.Examples
using DataFrames
using RCall
using Test
using Random: seed!

import PartiallyObservedMarkovProcesses as POMP

println("- Brownian motion tests")

@testset verbose=true "Brownian motion" begin

    seed!(1222664189)

    P = brownian_motion(
        times=[0:1:1000;5000:1:6000],
        x₀=[0,0],
        σ=[1.0 0.0;1.0 2.0],
        τ=[5.0 0; 1.0 1.0],
    );
    @test isa(P,POMP.PompObject)
    @time simulate(P,nsim=1000);
    @time simulate(P,nsim=1000);
    simulate_array(P,nsim=1000);
    @time simulate_array(P,nsim=1000);
    @time simulate_array(P,nsim=1000);

    d = transform(
        melt(P),
        :x => (x -> stack(x)[1,:]) => :x1,
        :x =>(x->stack(x)[2,:]) => :x2,
        :y => (x -> stack(x)[1,:]) => :y1,
        :y =>(x->stack(x)[2,:]) => :y2,
    );

    R"""
options(tidyverse.quiet=TRUE)
library(tidyverse)
$d |>
select(time,x1,x2,y1,y2) |>
filter(time<2000) |>
pivot_longer(-time) |>
ggplot(aes(x=time,y=value,group=name,color=name))+
geom_line()+
labs(color="")+
theme_bw()
"""

    R"""ggsave(filename="brown-01.png",width=7,height=4)"""

    R"""
options(tidyverse.quiet=TRUE)
library(tidyverse)
$d |>
select(time,x1,x2,y1,y2) |>
ggplot(aes(x=x1,y=x2))+
geom_point(alpha=0.3)+
geom_path(alpha=0.5)+
theme_bw()
"""

    R"""ggsave(filename="brown-02.png",width=7,height=7)"""

    @test_throws "different sizes" brownian_motion(
        times=[0:1:1000;5000:1:6000],
        x₀=[0,0],
        σ=[1.0 0.0;1.0 2.0],
        τ=[5.0 0 0; 1.0 1.0 0],
    )
    @test_throws "should be square" brownian_motion(
        times=[0:1:1000;5000:1:6000],
        x₀=[0,0],
        σ=[1.0 0.0 0;1.0 2.0 0],
        τ=[5.0 0 0; 1.0 1.0 0],
    )
    @test_throws "size mismatch" brownian_motion(
        times=[0:1:1000;5000:1:6000],
        x₀=[0,0,0],
        σ=[1.0 0.0;1.0 2.0],
        τ=[5.0 0; 1.0 1.0],
    )
    @test_throws "should be square" brownian_motion(
        times=[0:1:1000;5000:1:6000],
        x₀=[0,0,0,0],
        σ=[1.0 0.0 0;1.0 2.0 0],
        τ=[5.0 0 0; 1.0 1.0 0],
    )

end
