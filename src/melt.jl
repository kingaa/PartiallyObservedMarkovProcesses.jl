export melt

melt(x::Nothing) = DataFrame()

melt(x::NamedTuple; margins...) = DataFrame(;merge(x,margins)...)

melt(x::Array{<:NamedTuple}; margins...) =
    hcat(allcombinations(DataFrame;margins...),DataFrame(vec(x)))

melt(x::AbstractPompObject; id...) =
    hcat(melt(obs(x);id...,time=times(x)),melt(states(x)))

melt(
    x::AbstractArray{<:AbstractPompObject,N},
    names::Vararg{Symbol,N},
) where N = let
    c = NamedTuple{Tuple(names)}.(Tuple.(CartesianIndices(axes(x))))
    reduce(vcat,[melt(x[i];c[i]...) for i ∈ eachindex(x)])
end

melt(x::AbstractArray{<:AbstractPompObject,N}) where N =
    melt(x,Tuple(Symbol("id",i) for i ∈ 1:N)...)

melt(x::PfilterdPompObject) =
    hcat(
        melt(pomp(x)),
        DataFrame(
            ess=x.eff_sample_size,
            cond_logLik=x.cond_logLik
        )
    )
