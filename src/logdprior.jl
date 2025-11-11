"""
    logdprior(object; params=coef(object))

`logdprior` is the workhorse for the evaluator of the log prior density.
"""
logdprior(
    object::AbstractPompObject;
    params::Union{P,AbstractVector{P}} = coef(object),
) where {P<:NamedTuple} = begin
    try
        params = val_array(params)
        ell = Array{LogLik}(undef,size(params))
        logdprior!(object,ell;params=params)
        ell
    catch e
        if hasproperty(e,:msg)
            error("in `logdprior`: " * e.msg)
        else
            throw(e)            # COV_EXCL_LINE
        end
    end
end

"""
    logdprior!(object, ell; params=coef(object))

`logdprior!` is the in-place version of the `logdprior` workhorse.
"""
logdprior!(
    object::AbstractPompObject,
    ell::AbstractArray{W,1};
    params::Union{P,AbstractVector{P}} = coef(object),
) where {W<:Real,P<:NamedTuple} = begin
    try
        params = val_array(params)
        logdprior_internal!(ell,pomp(object).logdprior,params)
    catch e
        if isa(e,UndefKeywordError)
            error("in `logdprior!`: parameter " * string(e.var) * " undefined.")
        elseif isa(e,MethodError)
            error("in `logdprior!`: no matching method for args " * string(e.args[1]))
        elseif hasproperty(e,:msg)
            error("in `logdprior!`: " * e.msg)
        else
            throw(e)            # COV_EXCL_LINE
        end
    end
end

logdprior_internal!(          # COV_EXCL_LINE
    ell::AbstractArray{W,1},
    f::Nothing,
    _...,
) where {W<:Real} = begin
    for i ∈ eachindex(ell)
        @inbounds ell[i] = W(0)
    end
end

logdprior_internal!(
    ell::AbstractArray{W,1},
    f::Function,
    params::AbstractVector{P},
) where {W<:Real,P<:NamedTuple} = begin
    for j ∈ eachindex(params)
        @inbounds ell[j] = f(;params[j]...)::W
    end
end
