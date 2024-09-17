export simulate

mutable struct SimPompObject
    data::Union{Vector{NamedTuple},Nothing}
    t0::Real
    time::Vector{Real}
    params::Union{NamedTuple,Nothing}
    rinit::Union{Function,Nothing}
    rmeasure::Union{Function,Nothing}
    rprocess::Union{Function,Nothing}
    state::Array{<:NamedTuple,N} where N
end

convert(::Type{PompObject},object::SimPompObject) = PompObject(
    object.data,
    object.t0,
    object.time,
    object.params,
    object.rinit,
    object.rmeasure,
    object.rprocess
)

"""
    simulate(object; args...)

`simulate` simulates the POMP.
`args...` can be used to modify or unset fields.
"""
simulate(
    object::PompObject;
    args...
) = begin
    try
        object = pomp(object;args...)
        x0 = rinit(object)
        x = rprocess(object,x0=x0)
        y = rmeasure(object,x=x)
        object.data = reshape(y,length(y))
        SimPompObject(
            object.data,
            object.t0,
            object.time,
            object.params,
            object.rinit,
            object.rmeasure,
            object.rprocess,
            x
        )
    catch e
        if isa(e,UndefKeywordError)
            error("in `simulate`: parameter " * string(e.var) * " undefined.")
        elseif hasproperty(e,:msg)
            error("in `simulate`: " * e.msg)
        else
            throw(e)
        end
    end
end
