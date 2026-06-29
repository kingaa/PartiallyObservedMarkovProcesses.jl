module PartiallyObservedMarkovProcesses

import PartiallyObservedMarkovProcesses as POMP
using Reexport: @reexport

## Time is the (abstract) type for times.
const Time = Union{Int64,Float64}
const RealTime = Float64

## LogLik is the type for log likelihoods, etc.
## We may at some point need more precision....
const LogLik = Float64

export pomp, AbstractPompObject, PompPlugin, paramsymbs
include("pomp.jl")

export times, timezero, obs, init_state, states, coef
include("helpers.jl")

@reexport using Random
@reexport using Serialization: serialize, deserialize
export @bake, @freeze
include("bake.jl")

export rinit, rinit!
include("rinit.jl")

export rprocess, rprocess!
include("rprocess.jl")

export euler, discrete_time, onestep
include("plugins.jl")

export vectorfield
include("flow.jl")

export rmeasure
include("rmeasure.jl")

export logdmeasure, logdmeasure!
include("logdmeasure.jl")

export rprior
include("rprior.jl")

export logdprior, logdprior!
include("logdprior.jl")

export simulate, simulate_array
include("simulate.jl")

export pfilter, logLik, eff_sample_size, cond_logLik
include("pfilter.jl")

export traj_match_objfun
include("trajmatch.jl")

export melt
include("melt.jl")

export logmeanexp
include("logmeanexp.jl")

include("reshape.jl")

export EulerMultinomial, rand, rand!
include("eulermultinomial.jl")

# submodule Examples
include("examples/Examples.jl")

end # module
