export brownian_motion

brownian_motion = function(
    ;times::AbstractVector{<:Real},
    t₀::Real = 0,
    x₀::AbstractVector{<:Real},
    σ::AbstractMatrix{<:Real},
    τ::AbstractMatrix{<:Real},
    )
    if size(σ)!=size(τ)
        error("in `brownian_motion`: matrices σ and τ are of different sizes.")
    end
    if size(σ,1)!=size(σ,2)
        error("in `brownian_motion`: matrices σ and τ should be square.")
    end
    if size(σ,1)!=length(x₀)
        error("in `brownian_motion`: size mismatch between σ and x₀.")
    end
    simulate(
        t0=Float64(t₀),
        times=collect(Float64,times),
        params=(
            σ=collect(Float64,σ),
            τ=collect(Float64,τ),
            x₀=collect(Float64,x₀),
        ),
        rinit=function(;x₀,_...)
            (x=x₀,)
        end,
        rprocess=onestep(
            function(;dt,x,σ,_...)
                Δx = sqrt(dt)*σ*randn(Float64,length(x))
                (x=x+Δx,)
            end
        ),
        rmeasure=function(;x,τ,_...)
            Δy = τ*randn(Float64,length(x))
            (y=x+Δy,)
        end,
        logdmeasure=function(;y,x,τ,_...)
            logpdf(MvNormal(x,transpose(τ)*τ),y)
        end,
        nsim=1
    )[1]
end
