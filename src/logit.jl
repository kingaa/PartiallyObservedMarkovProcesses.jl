"""
    logit(p)

Logistic function: `logit(p) = log(p/(1-p))`.
See also [`expit`](@ref `expit`).
"""
@generated logit(p) = :(log(p/(1-p)))

"""
    expit(x)

Inverse logistic function: `expit(x) = 1/(1+exp(-x))`.
See also [`logit`](@ref `logit`).
"""
@generated expit(x) = :(1/(1+exp(-x)))

"""
    barycentric([type], x, n = 1)

Project Euclidean coordinates `x` onto a simplex.  Specifically, if `y
= barycentric(x,n)`, then `∑ y = n`.  Optionally, the results are
rounded to the nearest representable `type`.
"""
@generated barycentric(
    x, n::Real = 1,
) = begin
    if x <: NamedTuple
        names = fieldnames(x)
        quote
            (;zip($names, barycentric(values(x), n))...)
        end
    else
        quote
            m = Float64(n)/sum(Float64.(x))
            m.*x
        end
    end
end

@generated barycentric(
    type::Type{D}, x, n::Real = 1,
) where {D <: Real} = begin
    if x <: NamedTuple
        names = fieldnames(x)
        quote
            (;zip($names, barycentric(type, values(x), n))...)
        end
    else
        quote
            m = Float64(n)/sum(Float64.(x))
            round.(type, m.*x)
        end
    end
end
