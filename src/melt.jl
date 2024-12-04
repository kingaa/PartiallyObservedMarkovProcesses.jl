export melt

import DataFrames: DataFrame

melt(x::Nothing) = DataFrame()

melt(x::NamedTuple; id...) = DataFrame(;merge(NamedTuple(id),x)...)

melt(x::AbstractPompObject; id...) = let
    timevar = pomp(x).timevar
    t = NamedTuple{(timevar,)}.(times(x))
    n = length(t)
    s = isnothing(states(x)) ? fill((;),n) : states(x)
    d = isnothing(obs(x)) ? fill((;),n) : obs(x)
    DataFrame(
        map(
            (t,d,s) -> merge((;id...),t,d,s),
            t,d,s
        )
    )
end

melt(x::PfilterdPompObject; id...) =
    hcat(
        melt(pomp(x);id...),
        DataFrame(
            ess=x.eff_sample_size,
            cond_logLik=x.cond_logLik
        )
    )

melt(
    x::AbstractArray{T,N},
    names::Vararg{Symbol,N},
) where {T,N} = let
    c = NamedTuple{Tuple(names)}.(Tuple.(CartesianIndices(axes(x))))
    reduce(vcat,[melt(x[i];c[i]...) for i ∈ eachindex(x)])
end

melt(x::AbstractArray{T,N}) where {T,N} =
    melt(x,(Symbol("id",i) for i ∈ 1:N)...)
