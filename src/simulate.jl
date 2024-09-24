export simulate, simulate!, states

"""
    simulate(object; args...)

`simulate` simulates the POMP.
`args...` can be used to modify or unset fields.
"""
simulate(
    object::AbstractPompObject;
    args...,
) = begin
    try
        object = pomp(object;args...)
        x0 = rinit(object)
        x = rprocess(object,x0=x0)
        y = rmeasure(object,x=x)
        statezero!(object,x0)
        states!(object,x)
        obs!(object,y)
        object
    catch e
        if hasproperty(e,:msg)
            error("in `simulate`: " * e.msg)
        else
            throw(e)            # COV_EXCL_LINE
        end
    end
end

"""
    simulate!(object; args...)

`simulate!` simulates in place.
`args...` can be used to modify or unset fields.
"""
simulate!(
    object::AbstractPompObject;
    args...,
) = begin
    try
        pomp!(object;args...)
        x0 = rinit(object)
        x = rprocess(object,x0=x0)
        y = rmeasure(object,x=x)
        statezero!(object,x0)
        states!(object,x)
        obs!(object,y)
        nothing
    catch e
        if hasproperty(e,:msg)
            error("in `simulate!`: " * e.msg)
        else
            throw(e)            # COV_EXCL_LINE
        end
    end
end
