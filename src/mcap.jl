import DataFrames: DataFrame
import Distributions: Chisq, quantile
import LinearAlgebra: dot, qr, UpperTriangular, SingularException

# Local quadratic regression, as R's `loess(degree = 2, surface = "direct")`
# with one predictor.
struct Loess
    x::Vector{Float64}
    y::Vector{Float64}
    span::Float64
end

loess(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}; span::Real = 0.75) = begin
    @assert length(x) == length(y) "`x` and `y` must have the same length"
    @assert length(x) ≥ 3 "at least three observations are needed"
    @assert span > 0 "`span` must be positive"
    @assert all(isfinite,x) && all(isfinite,y) "`x` and `y` must be finite"
    # a window of q points gives the farthest weight zero, so at most q-1
    # carry weight, and a local quadratic needs 3 distinct values
    @assert span > 1 || floor(Int,span*length(x)) ≥ 4 "`span*n` must be at least 4"
    Loess(collect(Float64,x),collect(Float64,y),Float64(span))
end

(fit::Loess)(x0::Real) = begin
    x, y, span = fit.x, fit.y, fit.span
    n = length(x)
    d = abs.(x .- x0)
    h = if span ≤ 1
        partialsort(d,floor(Int,span*n))
    else
        # as R: the largest distance times sqrt(span)
        maximum(d)*sqrt(span)
    end
    @assert h > 0 "zero bandwidth at $x0: too few distinct `x` values"
    w = @. ifelse(d < h, (1-(d/h)^3)^3, 0.0)
    # A local quadratic needs at least three distinct x values with
    # positive weight; with fewer it is not determined by the data (R then
    # falls back on a pseudoinverse and warns), so refuse.
    @assert length(unique(x[w .> 0])) ≥ 3 "`span` = $(span) is too small: near $(x0) the smoother's window holds fewer than three distinct parameter values; increase `span`"
    # Weighted quadratic in (x-x0)/h. The intercept is the prediction;
    # scaling by h keeps the fit insensitive to the units of x.
    u = (x .- x0)./h
    X = hcat(ones(n),u,u.^2)
    sw = sqrt.(w)
    β = (sw .* X) \ (sw .* y)
    β[1]
end

"""
    MCAP

Created by a call to [`mcap`](@ref), this struct holds the Monte Carlo
adjusted profile.
"""
struct MCAP
    "profile log likelihoods"
    logLik::Vector{Float64}
    "profiled parameter values"
    parameter::Vector{Float64}
    "confidence level"
    level::Float64
    "smoothing span"
    span::Float64
    "evaluation grid, with the smoothed and quadratic fits"
    fit::DataFrame
    "grid maximizer of the smoothed profile"
    mle::Float64
    "maximizer of the quadratic fit (NaN if it is not concave)"
    quadratic_max::Float64
    "Monte Carlo adjusted confidence interval"
    ci::Tuple{Float64,Float64}
    "log-likelihood cutoff defining the interval"
    delta::Float64
    "statistical standard error"
    se_stat::Float64
    "Monte Carlo standard error"
    se_mc::Float64
    "total standard error"
    se::Float64
    "quadratic fit c + b*parameter - a*parameter^2"
    coefs::NamedTuple{(:c,:a,:b),NTuple{3,Float64}}
end

Base.show(io::IO, m::MCAP) = begin
    print(io,"MCAP: mle=$(m.mle), ci=$(m.ci), se=$(m.se) (stat $(m.se_stat), mc $(m.se_mc)), level=$(m.level)")
end

"""
    mcap(loglik, parameter; level = 0.95, span = 0.75, Ngrid = 1000)

Monte Carlo adjusted profile (Ionides et al. 2017), as R `pomp`'s
`mcap`.  `loglik` and `parameter` are typically the `loglik` column
and the profiled column returned by [`profile`](@ref); `span` must lie
in (0,1].  Returns an [`MCAP`](@ref).  The smoother corresponds to R's
`loess` with `surface = "direct"`, which differs slightly from R's
default.  Unlike R's, if the quadratic fit is not concave, the
standard errors, cutoff, and interval are `NaN`.  Each local fit needs
at least three distinct parameter values with positive weight: in the
quadratic window, everything derived from the quadratic is otherwise
`NaN` (R's `NA`), and in the smoother's windows the call is refused, as
the fit is not determined by the data.  In practice a profile needs at
least five distinct parameter values.  Rows whose log likelihood is
-∞ must be removed first; removing them biases the smoother upward at
those values.
"""
mcap(
    loglik::AbstractVector{<:Real},
    parameter::AbstractVector{<:Real};
    level::Real = 0.95,
    span::Real = 0.75,
    Ngrid::Integer = 1000,
) = begin
    ll = collect(Float64,loglik)
    par = collect(Float64,parameter)
    n = length(ll)
    @assert n == length(par) "`loglik` and `parameter` must have the same length"
    @assert all(isfinite,ll) "`loglik` must be finite"
    @assert all(isfinite,par) "`parameter` must be finite"
    @assert 0 < level < 1 "`level` must lie in (0,1)"
    @assert Ngrid ≥ 2 "`Ngrid` must be at least 2"
    @assert 0 < span ≤ 1 "`span` must lie in (0,1] for `mcap`"
    @assert trunc(Int,span*n) ≥ 1 "`span*length(parameter)` must be at least 1"
    smooth_fit = loess(par,ll;span)
    grid = collect(range(minimum(par),maximum(par),length=Ngrid))
    smoothed = smooth_fit.(grid)
    smooth_arg_max = grid[argmax(smoothed)]
    dist = abs.(par .- smooth_arg_max)
    included = dist .< sort(dist)[trunc(Int,span*n)]
    @assert any(included) "no points fall inside the quadratic fitting window; increase `span`"
    maxdist = maximum(dist[included])
    w = zeros(Float64,n)
    w[included] .= (1 .- (dist[included]./maxdist).^3).^3
    # A quadratic needs at least three distinct parameter values with
    # positive weight (the farthest included point has weight zero).  With
    # fewer it is not determined by the data: as with R's NA, everything
    # derived from it is NaN.
    ndistinct = length(unique(par[w .> 0]))
    determined = ndistinct ≥ 3
    determined || @warn "`mcap`: only $ndistinct distinct parameter values carry weight in the quadratic fit, so it is not determined by the data; its coefficients, standard errors, and interval are NaN. Increase `span` or add points."
    # weighted least squares, ll ~ c - a*z^2 + b*z, in the standardized
    # parameter z = (parameter-m)/s; results are converted back below
    m = (maximum(par)+minimum(par))/2
    s = (maximum(par)-minimum(par))/2
    z = (par .- m)./s
    X = hcat(ones(n),-z.^2,z)
    sw = sqrt.(w)
    F = qr(sw .* X)
    β = determined ? F \ (sw .* ll) : fill(NaN,3)
    c, a, b = β
    r = ll .- X*β
    nnz = count(>(0),w)
    # With no residual degrees of freedom or a singular design,
    # the standard errors are undefined.
    (nnz > 3 || !determined) || @warn "`mcap`: only $nnz points carry weight in the quadratic fit; increase `span` or add points."
    var_a,var_b,cov_ab = if determined && nnz > 3
        try
            σ2 = dot(w,r.^2)/(nnz-3)
            Ri = inv(UpperTriangular(F.R))
            V = σ2 .* (Ri*Ri')
            V[2,2],V[3,3],V[2,3]
        catch e
            e isa Union{ArgumentError,SingularException} || rethrow()
            @warn "`mcap`: the weighted design is singular; standard errors are not defined."
            NaN,NaN,NaN
        end
    else
        NaN,NaN,NaN
    end
    concave = determined && a > 0
    (concave || !determined) || @warn "`mcap`: the quadratic fit is not concave; standard errors and the interval are not defined."
    # standard errors in z units; `delta` is the same in any units.
    # All assume a concave fit, so otherwise they are NaN (R's are not)
    se_mc_squared = concave ? (1/(4*a*a))*(var_b - (2*b/a)*cov_ab + (b*b/a/a)*var_a) : NaN
    se_stat_squared = concave ? 1/2/a : NaN
    se_total_squared = se_mc_squared + se_stat_squared
    delta = concave ? quantile(Chisq(1),level)*(a*se_mc_squared + 0.5) : NaN
    logLik_diff = maximum(smoothed) .- smoothed
    inside = grid[logLik_diff .< delta]
    # the cutoff assumes a concave fit; otherwise there is no interval
    ci = (concave && !isempty(inside)) ? (minimum(inside),maximum(inside)) : (NaN,NaN)
    # as R: the square root of a negative variance is NaN
    nsqrt(v) = v ≥ 0 ? sqrt(v) : NaN
    zg = (grid .- m)./s
    quadratic = c .+ b.*zg .- a.*zg.^2
    MCAP(
        ll,par,Float64(level),Float64(span),
        DataFrame(parameter=grid,smoothed=smoothed,quadratic=quadratic),
        smooth_arg_max,concave ? m+s*b/(2*a) : NaN,ci,delta,
        s*nsqrt(se_stat_squared),s*nsqrt(se_mc_squared),s*nsqrt(se_total_squared),
        (c=c-b*m/s-a*m^2/s^2, a=a/s^2, b=b/s+2*a*m/s^2),
    )
end

mcap(_...) = error("Incorrect call to `mcap`.")
