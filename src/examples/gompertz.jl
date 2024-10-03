export gompertz

"""
    gompertz()

`gompertz` is a *PompObject* containing *Parus major* data and a simple Gompertz population model.

## Parameters
- r: the growth rate
- K: the equilibrium population density
- x₀: the initial population density
- σₚ: process noise s.d.
- σₘ: measurement noise s.d.
"""
gompertz = function()
    pomp(
        parus_data,
        t0=1960,
        times=:year,
        rinit = function (;x₀,_...)
            (;x=x₀,)
        end,
        rprocess = function (;t,x,σₚ,r,K,_...)
            s = exp(-r)
            d = LogNormal(s*log(x)+(1-s)*log(K),σₚ)
            (;t=t+1,x=rand(d),)
        end,
        rmeasure = function (;x,σₘ,_...)
            d = LogNormal(log(x),σₘ)
            (;pop=rand(d),)
        end,
        logdmeasure = function (;pop,x,σₘ,_...)
            logpdf(LogNormal(log(x),σₘ),pop)
        end
    )
end
