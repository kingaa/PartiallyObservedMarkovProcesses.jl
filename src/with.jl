using DataFrames: DataFrame

"""
    @with vals expr

Evaluates `expr` in a scope within which the values determined by the
named tuple `vals` are defined.
"""
macro with(vals, expr)
    esc(
        quote
            let vals = $vals, expr = $(QuoteNode(expr))
                local vlist
                vlist = PartiallyObservedMarkovProcesses.withlist(vals)
                eval(Expr(:let,vlist,expr))
            end
        end
    )
end

withlist(v::Union{NamedTuple,AbstractDict}) = begin
    Expr(:block,map((x,y) -> :($x = $y), keys(v), values(v))...)
end

withlist(v::DataFrame) = begin
    Expr(
        :block,
        map((x,y) -> :($x = $y),
            propertynames(v),
            collect(eachcol(v))
            )...,
    )
end

withlist(v) = error("@with will not work with `vals` of type $(typeof(v)).")

