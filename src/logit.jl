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
    barycentric(x)

Project Euclidean coordinates onto the unit simplex.
"""
@generated barycentric(x, n = 1) = quote
    m = n/sum(x)
    n.*x
end

barycentric(x::X, n = 1) where {X <: NamedTuple} = begin
    (;zip(fieldnames(X),barycentric(values(x),n))...)
end
