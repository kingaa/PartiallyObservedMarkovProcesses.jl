using Distributions: Normal, LogNormal, LogitNormal

macro normal(param, sd)
    draw = :($param = rand(Normal($param, scale * $sd)))
    (;param=param,draw=draw,)
end

macro lognormal(param, sd)
    draw = :($param = $param * rand(LogNormal(0.0, scale * $sd)))
    (;param=param,draw=draw,)
end

macro logitnormal(param, sd)
    draw = :($param = rand(LogitNormal(logit($param), scale * $sd)))
    (;param=param,draw=draw,)
end

macro logbarynormal(params, sd::Real)
    draws = map(params.args) do param
        :($param = $param * (rand(LogNormal(0, scale * $sd))))
    end
    names = Expr(:tuple,Expr(:parameters,params.args...))
    body = quote
        $names = begin
            $(draws...)
            barycentric($names)
        end
        $names
    end
    (;param=params.args,draw=body,)
end

macro logbarynormal(params, sd)
    draws = broadcast(params.args, sd.args) do param,s
        :($param = $param * rand(LogNormal(0, scale * $s)))
    end
    names = Expr(:tuple,Expr(:parameters,params.args...))
    body = quote
        $names = begin
            $(draws...)
            barycentric($names)
        end
        $names
    end
    (;param=params.args,draw=body,)
end

"""
    @ivp(expr, lag = 0)

Construct an initial-value parameter perturbation. The perturbation will be applied only at the specified lag. Thus, by default, the perturbation is applied at the zero-time (`lag = 0`) and if `lag = k`, then it will be applied after the `k`-th observation (but never after the final observation).
**This facility is experimental: the interface may change without warning.**
"""
macro ivp(expr, lag = 0)
    expr = eval(expr)
    param = flatten(expr.param)
    names = if length(param) > 1
        Expr(:tuple,Expr(:parameters,param...))
    else
        param[1]
    end
    draw = expr.draw
    body = quote
        $names = if lag == $lag
            $draw
        else
            $names
        end
    end
    (;param=param,draw=body,)
end

flatten() = []
flatten(a::Vector, b...) = [flatten(a...)..., flatten(b...)...]
flatten(a::Symbol, b...) = [a, flatten(b...)...]

"""
    @perturbn

Constructs a function suitable for use as the `perturbations` argument in [`mif`](@ref `mif`). For example, the following corresponds to log-normal perturbations of parameter `β`, logit-normal perturbations of parameter `p`, and log-barycentric-normal perturbations of `S₀,I₀,R₀`.  The latter is treated as an initial-value parameter.  **This facility is experimental: the interface may change without warning.**
```
    @perturbn(
        @lognormal(β,0.02),
        @logitnormal(p,0.02),
        @ivp(@logbarynormal((S₀,I₀,R₀),0.1))
    )
```
"""
macro perturbn(pieces...)
    pv = eval.([pieces...])
    args = flatten(map(x->x.param,pv))
    draws = map(pv) do piece
        quote
            $(piece.draw)
        end
    end
    names = Expr(:tuple,Expr(:parameters,args...))
    @eval function (scale, lag; $(args...), _...)
        $(draws...)
        $names
    end
end
