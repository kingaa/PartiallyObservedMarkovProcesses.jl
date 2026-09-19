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
