time_vector(t::Vector{Real}) = t
time_vector(t::Real) = [t]
time_vector(t,_...) = throw(ArgumentError("`time_vector` is not defined for type '"*string(typeof(t))*"'."))

param_vector(p::Vector{<:NamedTuple}) = p
param_vector(p::NamedTuple) = [p]
param_vector(p,_...) = throw(ArgumentError("`param_vector` is not defined for type '"*string(typeof(p))*"'."))

state_array(x::Array{<:NamedTuple,3},m::Integer,n::Integer) =
    if (size(x,2) != m || size(x,3) != n)
        error("state-array size mismatch")
    else
        x
    end
state_array(x::NamedTuple,_...) = reshape([x],1,1,1)
state_array(x,_...) = throw(ArgumentError("`state_array` is not defined for type '"*string(typeof(x))*"'."))
