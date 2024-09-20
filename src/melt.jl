export melt

melt(x::Array{<:NamedTuple,N};margins...) where N = 
    hcat(allcombinations(DataFrame,margins...),DataFrame(x))

melt(x::AbstractPompObject) =
    melt(obs(x),time=times(x))

melt(x::AbstractSimPompObject) =
    hcat(melt(states(x),time=times(x)),DataFrame(obs(x)))
