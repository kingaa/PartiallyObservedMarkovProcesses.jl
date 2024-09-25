export systematic_resample

systematic_resample(
    w::Array{Float64,N},
    Np::Int64,
) where N = begin
    w = cumsum(vec(w))
    n = length(w)
    @assert w[n] > 0 "in `systematic_resample`: sum of weights should be positive!"
    p = Array{Int64}(undef,Np)
    du = w[n]/Float64(Np);
    u = -du*Random.rand()
    i = 1
    for j ∈ 1:Np
        u += du
        while (u > w[i] && i < n)
            i += 1
        end
        p[j] = i
    end
    p                           # COV_EXCL_LINE
end
