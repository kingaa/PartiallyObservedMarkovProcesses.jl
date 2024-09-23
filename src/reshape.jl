vectorize(t::Vector{T}) where T = t

vectorize(t) = [t]

val_array(x::Union{Vector{<:NamedTuple},Array{<:NamedTuple,N}},dim...) where N = begin
    n = prod(dim)
    if mod(length(x),n) != 0
        error("size mismatch in `val_array`")
    else
        reshape(x,dim...,div(length(x),n))
    end
end

val_array(x::NamedTuple,_...) = reshape([x],1,1,1)

val_array(x,_...) = throw(ArgumentError("`val_array` is not defined for type '"*string(typeof(x))*"'."))
