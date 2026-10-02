using PartiallyObservedMarkovProcesses
using PartiallyObservedMarkovProcesses.Examples
using Random
using Statistics: mean, std
using Test

@info h1("Gompertz model: exact Kalman likelihood check")

# The Gompertz state-space model reduces to a linear-Gaussian state-space
# model on the log scale, so its likelihood can be computed exactly (up to
# a change-of-variables constant) by the scalar Kalman filter, with no
# Monte Carlo error. This gives an exact target against which the
# particle filter's log-likelihood estimate can be checked.
#
# Writing Z_t = log(X_t), the process model
#   X_t = X_{t-1}^S K^(1-S) eps_t,   S = exp(-r),   eps_t ~ LogNormal(0,σₚ),
# becomes the linear-Gaussian autoregression
#   Z_t = (1-S) log(K) + S Z_{t-1} + w_t,   w_t ~ Normal(0,σₚ²),
# and the measurement model
#   pop_t ~ LogNormal(Z_t,σₘ)
# becomes, on the log scale,
#   log(pop_t) = Z_t + v_t,   v_t ~ Normal(0,σₘ²).
# `rinit` sets X_0 to the parameter X0 with no added noise, and the first
# observation time in the built-in data coincides with t0, so no process
# step intervenes before the first observation: the filtering
# distribution of Z at the first observation time is a point mass at
# log(X0), i.e. a Normal with zero variance. Every subsequent step
# applies one Gompertz recursion (dt=1) before the next observation.
# The `discrete_time` step function itself ignores its `dt` argument
# (it uses S=exp(-r) and passes σₚ bare to the LogNormal, with no
# sqrt(dt) rescaling), so σₚ is exactly the per-step log-scale standard
# deviation whenever, as here, `discrete_time` is called with dt=1 over
# unit-spaced annual observation times.
#
# The following is the textbook scalar Kalman filter recursion (predict,
# then update on the Gaussian innovation), returning the exact
# conditional log likelihood of log(pop_1),...,log(pop_N) under this
# linear-Gaussian state-space model.
kalman_loglik_logy = function (logy,r,K,σₚ,σₘ,X0)
    S = exp(-r)
    a = (1-S)*log(K)
    Q = σₚ^2
    R = σₘ^2
    m = log(X0)   # filtering mean of Z_1, before assimilating obs 1
    v = 0.0        # filtering variance of Z_1: point mass, since X_1=X0
    ll = 0.0
    for t ∈ eachindex(logy)
        if t > 1
            # predict: one Gompertz step of the linear-Gaussian recursion
            m = a + S*m
            v = S^2*v + Q
        end
        # innovation and its (exactly Gaussian) conditional log likelihood
        fvar = v + R
        innov = logy[t] - m
        ll += -0.5*(log(2π*fvar) + innov^2/fvar)
        # update: filtering distribution of Z_t given obs up to time t
        gain = v/fvar
        m = m + gain*innov
        v = v - gain^2*fvar
    end
    ll
end

@testset verbose=true "Gompertz: exact Kalman likelihood check" begin

    P = gompertz()
    p1 = (r=4.5,K=210.0,σₚ=0.7,σₘ=0.1,X0=150.0)
    y = [o.pop for o ∈ obs(P)]
    logy = log.(Float64.(y))

    # exact conditional log likelihood of \log(pop_{1:N}) under the
    # linear-Gaussian reduction of the Gompertz state-space model
    ll_logy = kalman_loglik_logy(logy,p1.r,p1.K,p1.σₚ,p1.σₘ,p1.X0)
    @test isfinite(ll_logy)

    # `logdmeasure` in the Gompertz model evaluates the LogNormal density
    # of pop_t itself, not of log(pop_t). Since pop_t=exp(log(pop_t)),
    # the two densities are related by the Jacobian of this change of
    # variables, f_pop(y) = f_logpop(log y)/y, so
    #   loglik(pop_1:N) = loglik(log(pop_1:N)) - sum_t log(pop_t).
    jacobian = -sum(logy)
    ll_y = ll_logy + jacobian
    @test isfinite(ll_y)

    # the particle filter (classical bootstrap filter: trigger=1,
    # target=0, the current defaults) estimates loglik(pop_1:N),
    # exactly this quantity, up to Monte Carlo error that vanishes as
    # Np grows.
    Random.seed!(20260907)
    Nps = [200,1000,5000]
    nreps = 30
    lme = [
        logmeanexp([logLik(pfilter(P,Np=Np,params=p1)) for _ ∈ 1:nreps],se=true)
        for Np ∈ Nps
    ]
    gap = [x.est - ll_logy for x ∈ lme]
    se = [x.se for x ∈ lme]

    # primary check: at every Np, the estimate agrees with the exact
    # log likelihood (the Kalman value plus the Jacobian) within Monte
    # Carlo error
    for (g,s) ∈ zip(gap,se)
        @test abs(g-jacobian) < 6*s
    end

    # supplementary check: the estimates at the different Np agree with
    # one another, within the Monte Carlo error at the smallest Np
    tol_spread = 6*maximum(se)
    @test maximum(gap)-minimum(gap) < tol_spread

end

@testset verbose=true "Gompertz: slice against the exact likelihood" begin

    # A likelihood slice holds every other parameter fixed, so each row
    # is a fixed-parameter likelihood and can be checked, point by
    # point, against the exact Kalman value plus the Jacobian constant.
    P = gompertz()
    p1 = (r=4.5,K=210.0,σₚ=0.7,σₘ=0.1,X0=150.0)
    y = [o.pop for o ∈ obs(P)]
    logy = log.(Float64.(y))
    jacobian = -sum(logy)
    d = slice_design(p1; K=[150.0,180.0,210.0,240.0,270.0], σₘ=[0.05,0.2])
    Random.seed!(20260922)
    s = slice(P,d;Np=2000,nreps=10)
    exact = [
        kalman_loglik_logy(logy,r.r,r.K,r.σₚ,r.σₘ,r.X0)+jacobian
        for r ∈ eachrow(s)
    ]
    z = (s.loglik .- exact)./s.se
    @info "slice vs exact, standardized gaps: $(round.(z,digits=2))"
    @test all(abs.(z) .< 5)
    # the slice is not a profile: the other parameters are untouched
    @test all(s.r .== p1.r) && all(s.X0 .== p1.X0)

end

@testset verbose=true "Gompertz: profile against the exact profile likelihood" begin

    # Profile over K with σₚ as the one nuisance parameter. The exact
    # profile, max over σₚ of the exact likelihood, is found by golden
    # section search on log σₚ. Three checks:
    # (i) each reported log likelihood is a fixed-parameter estimate at
    #     exactly the reported parameter vector;
    # (ii) no reported value exceeds the exact profile beyond Monte
    #      Carlo error, since the profile is a maximum;
    # (iii) mif has come close to the maximizing σₚ.
    P = gompertz()
    p1 = (r=4.5,K=210.0,σₚ=0.7,σₘ=0.1,X0=150.0)
    y = [o.pop for o ∈ obs(P)]
    logy = log.(Float64.(y))
    jacobian = -sum(logy)
    exact(θ) = kalman_loglik_logy(logy,θ.r,θ.K,θ.σₚ,θ.σₘ,θ.X0)+jacobian
    golden(f,a,b) = begin
        g = (sqrt(5)-1)/2
        for _ ∈ 1:80
            lo = b-g*(b-a); hi = a+g*(b-a)
            f(lo) > f(hi) ? (b = hi) : (a = lo)
        end
        (a+b)/2
    end
    exact_profile(K) = begin
        θ(s) = merge(p1,(K=K,σₚ=exp(s)))
        exact(θ(golden(s->exact(θ(s)),log(0.01),log(3.0))))
    end

    d = profile_design(
        K=[150.0,180.0,210.0,240.0,270.0];
        lower=(σₚ=0.2,),upper=(σₚ=1.5,),
        nprof=2,rng=MersenneTwister(7),
    )
    d.r .= p1.r; d.σₘ .= p1.σₘ; d.X0 .= p1.X0
    Random.seed!(20260922)
    pr = profile(
        P,d;
        Nmif=40,Np=1000,
        perturbations=@perturbn(σₚ ~ LogNormal(0.1)),
        cooling=geometric_cooling(0.5),
        nreps=10,Np_eval=2000,
    )
    at = [exact(NamedTuple(r[[:r,:K,:σₚ,:σₘ,:X0]])) for r ∈ eachrow(pr)]
    ep = exact_profile.(pr.K)
    z = (pr.loglik .- at)./pr.se
    @info "profile vs exact at reported point, standardized gaps: $(round.(z,digits=2))"
    @info "exact profile minus exact value at reported point: $(round.(ep .- at,digits=3))"
    # (i) evaluated at the reported point, K and fixed parameters exact
    @test all(abs.(z) .< 5)
    @test pr.K == d.K && all(pr.r .== p1.r) && all(pr.X0 .== p1.X0)
    # (ii) the profile is a maximum
    @test all(pr.loglik .< ep .+ 5 .* pr.se)
    # (iii) mif has converged to within one log-likelihood unit
    @test all(0 .≤ ep .- at .< 1.0)

end

@testset verbose=true "Gompertz: likelihood-scale checks against the exact likelihood for selected resampling settings" begin

    # Triggered and power-tempered resampling should change the variance
    # of the likelihood estimate, not its expectation: for each setting
    # tested, the mean of the likelihood estimate itself (not of its
    # logarithm) is compared with the exact likelihood. σₘ = 0.3 keeps the effective sample
    # size high enough that a trigger below 1 really does skip
    # resampling at some observations; this is checked too, since with
    # an informative measurement model every setting resamples at every
    # step and the test would be blind to triggering.
    P = gompertz()
    p1 = (r=4.5,K=210.0,σₚ=0.3,σₘ=0.3,X0=150.0)
    logy = log.(Float64.([o.pop for o ∈ obs(P)]))
    exact = kalman_loglik_logy(logy,p1.r,p1.K,p1.σₚ,p1.σₘ,p1.X0)-sum(logy)
    Random.seed!(20260923)
    R = 1000
    for (trigger,target) ∈ ((1.0,0.0),(0.5,0.0),(0.5,0.5),(0.2,0.75))
        pfs = [pfilter(P;Np=200,params=p1,trigger,target) for _ ∈ 1:R]
        frac = mean(mean(resampled(pf)[2:end]) for pf ∈ pfs)
        ratio = [exp(logLik(pf)-exact) for pf ∈ pfs]
        z = (mean(ratio)-1)/(std(ratio)/sqrt(R))
        @info "trigger=$trigger, target=$target: resampled at $(round(100*frac))% of steps; mean of L̂/L = $(round(mean(ratio),digits=4)), z = $(round(z,digits=2))"
        @test abs(z) < 5
        trigger < 1 && @test 0 < frac < 1
    end

end
