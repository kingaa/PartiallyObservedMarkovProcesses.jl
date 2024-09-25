export melt

melt(x::Array{<:NamedTuple,N};margins...) where N = 
    hcat(allcombinations(DataFrame,margins...),DataFrame(x))

melt(x::AbstractPompObject) = melt(obs(x),time=times(x))

melt(x::SimPompObject) = begin
    t = times(x)
    s = states(x)
    o = obs(x)
    hcat(
        melt(s,rep=axes(s,1),parset=axes(s,2),time=t),
        melt(o)
    )
end
