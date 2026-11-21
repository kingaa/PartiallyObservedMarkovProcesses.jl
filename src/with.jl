using DataFrames: DataFrame

"""
    @with vals expr

Evaluates `expr` in a scope within which the variables determined by
`vals` are defined. The latter can be a `NamedTuple`, `Dict`, or
`DataFrame`.
"""
macro with(vals, expr)
    esc(
        quote
            let vals = $vals, expr = $(QuoteNode(expr))
                local vlist
                vlist = PartiallyObservedMarkovProcesses.withlist(vals)
                eval(Expr(:let, vlist, expr))
            end
        end
    )
end

withlist(v::Union{NamedTuple,AbstractDict}) = begin
    vlist = map((x,y) -> :($x = $y), keys(v), values(v))
    Expr(:block, vlist...)
end

withlist(v::DataFrame) = begin
    vlist = map(
        (x,y) -> :($x = $y),
        propertynames(v),
        collect(eachcol(v))
    )
    Expr(:block, vlist...)
end

withlist(v) = error("@with will not work with `vals` of type $(typeof(v)).")
