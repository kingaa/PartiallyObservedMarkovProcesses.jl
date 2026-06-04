import Statistics: std

logmeanexp1(x) = begin
    xmax = maximum(x)
    xmax + log(sum(exp.(x .- xmax))) - log(length(x))
end

logmeanexp1(x, drop) = begin
    xv = @views vcat(x[begin:drop-1], x[drop+1:end])
    logmeanexp(xv)
end

ess1(x) = begin
    w = exp.(x .- maximum(x))
    sum(w)^2/sum(w.^2)
end

"""
    logmeanexp(x; se = false, ess = false)

Compute the log-mean-exp of `x`. Optionally, return a jack-knife estimate
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
        jk = [logmeanexp1(x,i) for i ∈ eachindex(x)]
        xse = (n-1)*std(jk)/sqrt(n)
        if ess
            (est=lme, se=xse, ess=ess1(x))
        else
            (est=lme, se=xse)
        end
    elseif ess
        (est=lme, ess=ess1(x))
    else
        lme
    end
end
