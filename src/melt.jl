export melt

melt(x::Array{<:NamedTuple,N};margins...) where N = 
    hcat(allcombinations(DataFrame,margins...),DataFrame(x))

melt(x::AbstractPompObject) = if isnothing(states(x))
    melt(obs(x),time=times(x))
else
    hcat(melt(obs(x),time=times(x)),DataFrame(states(x)))
end    
