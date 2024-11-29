export systematic_resample, systematic_resample!

systematic_resample!(
    p::AbstractVector{I},
    w::AbstractVector{F},
) where {I<:Integer,F<:Real} = let
    Np = length(p)
    w = cumsum(w)
    n = length(w)
    @assert w[n] >= 0 "in `systematic_resample`: sum of weights should be positive!"
    du = w[n]/F(Np);
    u = -du*rand(F)
    i = 1
    for j ∈ eachindex(p)
        u += du
        @inbounds while (u > w[i] && i < n)
            i += 1
        end
        @inbounds p[j] = i
    end
end

systematic_resample(
    Np::Integer,
    w::AbstractVector{<:Real},
    I::Type = Int64,
) = let
    p = Vector{I}(undef,Np)
    systematic_resample!(p,w)
    p
end
