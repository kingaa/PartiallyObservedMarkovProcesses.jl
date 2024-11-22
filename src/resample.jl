export systematic_resample, systematic_resample!

systematic_resample!(
    p::AbstractVector{Int64},
    w::AbstractArray{Float64},
) = let
    Np = length(p)
    w = cumsum(vec(w))
    n = length(w)
    @assert w[n] >= 0 "in `systematic_resample`: sum of weights should be positive!"
    du = w[n]/Float64(Np);
    u = -du*Random.rand()
    i = 1
    for j ∈ eachindex(p)
        u += du
        while (u > w[i] && i < n)
            i += 1
        end
        p[j] = i
    end
end

systematic_resample(
    Np::Int64,
    w::AbstractArray{Float64},
) = let
    p = Array{Int64}(undef,Np)
    systematic_resample!(p,w)
    p                           # COV_EXCL_LINE
end
