export gompertz

import Distributions: LogNormal, logpdf

@doc raw"""
    gompertz()

`gompertz` is a *PompObject* containing *Parus major* data and a simple Gompertz population model.
The population model has a single scalar state variable, ``X_t``, which obeys
```math
X_t = X_{t-1}^S\,K^{1-S}\,\varepsilon_t,
```
where ``S = e^{-r\delta{t}}`` and ``\varepsilon_t \sim \mathrm{LogNormal}(0,\sigma_p)``.
The time-step is one unit: ``\delta{t}=1``.
The data are assumed to be drawn from a log-normal distribution.
In particular,
```math
\mathrm{pop}_t \sim \mathrm{LogNormal}(\log{X_t},\sigma_m).
```

## Parameters
- r: the growth rate
- K: the equilibrium population density
- X₀: the initial population density
- σₚ: process noise s.d.
- σₘ: measurement noise s.d.
"""
gompertz = function()
    pomp(
        parus_data,
        t0=1960,
        times=:year,
        rinit = function (;X₀,_...)
            (;X=X₀,)
        end,
        rprocess = euler(
            function (;t,X,σₚ,r,K,_...)
                s = exp(-r)
                d = LogNormal(s*log(X)+(1-s)*log(K),σₚ)
                (;X=rand(d),)
            end,
            dt=1
        ),
        rmeasure = function (;X,σₘ,_...)
            d = LogNormal(log(X),σₘ)
            (;pop=rand(d),)
        end,
        logdmeasure = function (;pop,X,σₘ,_...)
            logpdf(LogNormal(log(X),σₘ),pop)
        end
    )
end
