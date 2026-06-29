import Random: default_rng
import Distributions: rand, Binomial

struct EulerMultinomial{N, I <: Integer, F <: Real, R <: AbstractVector{F}}
    size::I
    rate::R
    lambda::F
    dt::Time
    EulerMultinomial(
        size::I, rate::R, dt::Time,
    ) where {I <: Integer, F <: Real, R <: AbstractVector{F}} = begin
        @assert isfinite(size) && size ≥ zero(I) "invalid size = $(size)"
        @assert isfinite(dt) && dt ≥ zero(Time) "invalid dt = $(dt)"
        lambda::F = zero(F)
        for k ∈ eachindex(rate)
            @assert isfinite(rate[k]) && rate[k] ≥ zero(F) "invalid rate[$k] = $(rate[k])"
            lambda += rate[k]
        end
        new{length(rate),I,F,R}(size,rate,lambda,RealTime(dt))
    end
end

"""
    rand!(trans, d::EulerMultinomial)

Fill the vector `trans` with a random draw from the Euler-multinomial distribution defined by `d`.
"""
rand!(
    rng::AbstractRNG,
    trans::AbstractVector{I},
    d::EulerMultinomial{N,I,F,R},
) where {N,I,F,R} = begin
    @assert length(trans) == N "size mismatch"
    p::F = d.lambda
    if p > zero(F) && d.size > zero(I)
        size::I = rand(Binomial(d.size,1-exp(-p * d.dt)))
        k::Int = 1
        while k < length(trans)
            trans[k] =  (size > zero(I) && p > zero(F)) ?
                rand(Binomial(size,d.rate[k]/p)) : zero(I)
            size -= trans[k]
            p -= d.rate[k]
            k += 1
        end
        trans[k] = size;
    else
        trans .= zero(I)
    end
    nothing
end

rand!(
    trans::AbstractVector{I},
    d::EulerMultinomial{N,I,F,R},
) where {N,I,F,R} = rand!(default_rng(),trans,d)

rand(
    rng::AbstractRNG,
    d::EulerMultinomial{N,I,F,R},
) where {N,I,F,R} = begin
    trans = Array{I}(undef,N)
    rand!(rng,trans,d)
    trans
end
