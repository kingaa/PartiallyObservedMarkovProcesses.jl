val_array(x::Array{X,N}) where {X,N} = vec(x)

val_array(x::Array{X,N}, dim...) where {X,N} = begin
    q,r = divrem(length(x),prod(dim))
    if r != 0
        error("in `val_array`: size mismatch.")
    else
        reshape(x,q,dim...)
    end
end

val_array(x) = [x]

val_array(x, dim...) = begin
    reshape([x],dim...)
end

