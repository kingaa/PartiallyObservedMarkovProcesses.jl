val_array(x::Vector) = x

val_array(x::Array{X,N}) where {X,N} = vec(x)

val_array(x::Array{X,N}, dim::Integer...) where {X,N} = begin
    q,r = divrem(length(x),prod(dim))
    if r != 0
        error("in `val_array`: size mismatch.")
    else
        reshape(x,dim...,q)
    end
end

val_array(x) = [x]

val_array(x, dim::Integer...) = val_array([x],dim...)
