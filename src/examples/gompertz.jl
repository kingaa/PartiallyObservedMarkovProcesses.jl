export gompertz

"""
    gompertz(r=4.5,K=210.0,σₚ=0.7,σₘ=0.1,x₀=150.0)

`gompertz` is a *PompObject* containing *Parus major* data and a simple Gompertz population model.

## Parameters
- r: the growth rate
- K: the equilibrium population density
- x₀: the initial population density
- σₚ: process noise s.d.
- σₘ: measurement noise s.d.
"""
gompertz = function(r::Float64=4.5, K::Float64=210.0, σₚ::Float64=0.7,
                    σₘ::Float64=0.1,
                    x₀::Float64=150.0)
    pomp(
        parus_data,
        t0=1960,
        times=:year,
        params=(;r=r,K=K,σₚ=σₚ,σₘ=σₘ,x₀=x₀),
        rinit = function (;x₀,_...)
            (;x=x₀,)
        end,
        rmeasure = function (;x,σₘ,_...)
            d = LogNormal(log(x),σₘ)
            (;pop=rand(d),)
        end,
        rprocess = function (;t,x,σₚ,r,K,_...)
            s = exp(-r)
            d = LogNormal(s*log(x)+(1-s)*log(K),σₚ)
            (;t=t+1,x=rand(d),)
        end
    )
end
