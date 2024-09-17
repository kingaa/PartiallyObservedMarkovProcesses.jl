time_vector(t::Vector{Real}) = t
time_vector(t::Real) = [t]
time_vector(t,_...) = throw(ArgumentError("`time_vector` is not defined for type '"*string(typeof(t))*"'."))

val_array(x::Union{Vector{<:NamedTuple},Array{<:NamedTuple,N}},dim...) where N = begin
    n = prod(dim)
    if (mod(length(x),n)!=0)
        error("state-array size mismatch")
    else
        reshape(x,div(length(x),n),dim...)
    end
end
val_array(x::NamedTuple,_...) = reshape([x],1,1,1)
val_array(x,_...) = throw(ArgumentError("`val_array` is not defined for type '"*string(typeof(x))*"'."))