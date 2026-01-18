import SciMLBase
import DifferentialEquations: ODEProblem, solve

struct VectorfieldPlugin{F<:Function, G<:Function} <: PompPlugin
    integrator::F
    vf_evaluator::G
    statenames::Vector{Symbol}
end

argnames(m::Method) = Base.rest(Base.method_argnames(m),2)

"""
    vectorfield(vf, integration_alg; integrator_args...)

The function `vf` should return a vector containing the
components of the vectorfield *in the same order* in which they are passed
to `vf`.

The `integration_alg` argument specifies the integration algorithm to be used.  See **DifferentialEquations.jl** for more on the choices.  The keyword arguments in `integrate_args...` are passed to the `DifferentialEquations.solve()` function.
"""
vectorfield(
    vf::Function,
    integration_alg;
    integrator_args...,
) = begin
    m = methods(vf)
    @assert length(m)==1
    statenames = argnames(m[1])
    vf_eval!(du, u, p, t) = begin
        du[:] = vf(u...;p...,t=t)
        nothing
    end
    integrator!(
        x::AbstractVector{X},
        t::AbstractVector{<:RealTime},
        x0::X,
        params::P,
    ) where {X<:NamedTuple,P<:NamedTuple} = begin
        tspan = extrema(t)
        ic = [x0[statenames]...]
        prob = ODEProblem{true,SciMLBase.NoSpecialize}(
            vf_eval!,
            ic,tspan,params
        )
        sol = solve(prob,integration_alg;integrator_args...)
        x[:] = map(v->(;zip(statenames,v)...),sol(t))
    end
    VectorfieldPlugin(integrator!,vf_eval!,statenames)
end

vectorfield(_...) = error("Incorrect call to `vectorfield`.")

## vectorfield (deterministic, continuous-time)
## advance the state for each IC and parameter
## without accumulators
rproc_internal!(
    x::AbstractArray{X,3},
    plugin::VectorfieldPlugin,
    x0::AbstractArray{X,2},
    times::AbstractVector{T},
    t0::T,
    params::AbstractVector{P},
    accumvars::Nothing,
    userdata::U,
) where {T<:RealTime,X<:NamedTuple,P<:NamedTuple,U<:NamedTuple} = begin
    for j ∈ eachindex(params), k ∈ axes(x0,2)
        @inbounds plugin.integrator(@view(x[:,j,k]),times,x0[j,k],(;params[j]...,userdata...))
    end
end
