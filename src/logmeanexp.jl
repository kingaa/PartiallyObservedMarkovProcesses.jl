using Statistics: std
using DataFrames: Not

logmeanexp1(x::AbstractVector{W}) where {W <: Real} = begin
    xmax = maximum(x)
    if isfinite(xmax)
        xmax + log(sum(exp.(x .- xmax))) - log(length(x))
    else
        xmax
    end
end

logmeanexp1(x, drop) = begin
    logmeanexp1(@view(x[Not(drop)]))
end

ess1(x) = begin
    xmax = maximum(x)
    w = exp.(x .- xmax)
    sum(w)^2/sum(w.^2)
end

"""
    logmeanexp(x; se = false, ess = false)

Compute `log(mean(exp(x)))`. Optionally, return a jack-knife estimate
 of the standard error (`se = true`) and/or the effective sample size
 (`ess = true`).
"""
logmeanexp(
    x;
    se = false,
    ess = false,
) = begin
    lme = logmeanexp1(x)
    if se
        n = length(x)
        jk = (n-1)*std(logmeanexp1(x,i) for i ∈ eachindex(x))/sqrt(n)
        if ess
            (est=lme, se=jk, ess=ess1(x))
        else
            (est=lme, se=jk)
        end
    elseif ess
        (est=lme, ess=ess1(x))
    else
        lme
    end
end
