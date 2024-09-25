val_array(x::Array{T,N},dim...) where {T,N} = begin
    n = prod(dim)
    if mod(length(x),n) != 0
        error("size mismatch in `val_array`")
    else
        reshape(x,div(length(x),n),dim...)
    end
end

val_array(x,dim...) = begin
    val_array([x],dim...)
end
