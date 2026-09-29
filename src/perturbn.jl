using Distributions: Normal, LogNormal, LogitNormal

normalpert(param, sd) = begin
    @assert(
        param isa Symbol,
        "proper specification is `param ~ Normal(sd)`"
    )
    draw = :($param += rand(Normal(0, scale * $sd)))
    (;param,draw,)
end

lognormalpert(param, sd) = begin
    @assert(
        param isa Symbol,
        "proper specification is `param ~ LogNormal(sd)`"
    )
    draw = :($param *= rand(LogNormal(0, scale * $sd)))
    (;param,draw,)
end

logitnormalpert(param, sd) = begin
    @assert(
        param isa Symbol,
        "proper specification is `param ~ LogitNormal(sd)`"
    )
    draw = :($param = rand(LogitNormal(logit($param), scale * $sd)))
    (;param,draw,)
end

logbarynormalpert(params, sd) = begin
    @assert(
        params isa Expr && params.head == :tuple && length(params.args) > 1,
        "proper specification is `(p1,p2,...) ~ LogBaryNormal(sd)`"
    )
    draws = map(params.args) do p
        @assert p isa Symbol "invalid parameter `$p`"
        :($p = $p * rand(LogNormal(0, scale * $sd)))
    end
    names = Expr(:tuple, Expr(:parameters, params.args...))
    draw = quote
        $(draws...)
        barycentric($names)
    end
    (;param=params.args,draw,)
end

disallowed_perturbation_type_check(expr, allowable) = begin
    @assert(
        expr.head == :call && expr.args[1] in allowable,
        begin
            allowed = join(map(s->"`"*string(s)*"`", allowable), ", ")
            """Unrecognized perturbation specification `$expr`: Perturbation type should be one of: $allowed."""
        end
    )
end

ivppert(param, expr, lag = 0) = begin
    disallowed_perturbation_type_check(
        expr,
        [:Normal, :LogNormal, :LogitNormal, :LogBaryNormal]
    )
    draws = randpert(param, expr)
    param = flatten(draws.param)
    names = if length(param) > 1
        Expr(:tuple, Expr(:parameters, param...))
    else
        param[1]
    end
    draw = quote
        $names = if lag == $lag
            $(draws.draw)
        else
            $names
        end
    end
    (;param,draw,)
end

randpert(param, (expr...,)) = begin
    disallowed_perturbation_type_check(
        expr,
        [:ivp, :Normal, :LogNormal, :LogitNormal, :LogBaryNormal]
    )
    type, rest... = expr.args
    if type==:ivp
        ivppert(param, rest...)
    elseif type==:Normal
        normalpert(param, rest...)
    elseif type==:LogNormal
        lognormalpert(param, rest...)
    elseif type==:LogitNormal
        logitnormalpert(param, rest...)
    elseif type==:LogBaryNormal
        logbarynormalpert(param, rest...)
    end
end

flatten() = []
flatten(a::Vector, b...) = [flatten(a...)..., flatten(b...)...]
flatten(a::Symbol, b...) = [a, flatten(b...)...]

"""
    @perturbn

Constructs a function suitable for use as the `perturbations` argument in [`mif`](@ref `mif`) according to specifications given using a simple mini-language. For example, the following corresponds to log-normal perturbations of parameter `β`, logit-normal perturbations of parameter `p`, and log-barycentric-normal perturbations of `(S₀,I₀,R₀)`.  The latter is treated as an initial-value parameter.  **This facility is experimental: the interface may change without warning.**
```
    @perturbn(
        β ~ LogNormal(0.02),
        p ~ LogitNormal(0.02),
        (S₀,I₀,R₀) ~ ivp(LogBaryNormal(0.1))
    )
```
"""
macro perturbn(pieces...)
    pv = map([pieces...]) do piece
        @assert(
            piece isa Expr,
            "proper `@perturbn` specification is `param ~ type(sd)` or `param ~ ivp(type(sd),[lag])."
        )
        op, param, rest... = piece.args
        @assert(
            piece.head==:call && op==:~,
            "proper `@perturbn` specification is `param ~ type(sd)` or `param ~ ivp(type(sd), [lag])."
        )
        randpert(param, rest...)
    end
    args = flatten(map(x->x.param, pv))
    draws = map(pv) do piece
        quote
            $(piece.draw)
        end
    end
    names = Expr(:tuple, Expr(:parameters, args...))
    @eval function (scale, lag; $(args...), _...)
        $(draws...)
        $names
    end
end
