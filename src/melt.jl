export melt

melt(x::Array{<:NamedTuple,N};margins...) where N =
    hcat(allcombinations(DataFrame,margins...),DataFrame(x))

melt(x::AbstractPompObject) = melt(obs(x),time=times(x))

melt(x::SimPompObject) = begin
    s = states(x)
    hcat(
        melt(s,rep=axes(s,1),parset=axes(s,2),time=times(x)),
        melt(obs(x))
    )
end

melt(x::PfilterdPompObject) = begin
    hcat(
        melt(obs(x),time=times(x)),
        DataFrame(
            ess=x.eff_sample_size,
            cond_logLik=x.cond_logLik
        )
    )
end
