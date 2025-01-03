module PartiallyObservedMarkovProcesses

import PartiallyObservedMarkovProcesses as POMP
import DataFrames: DataFrame, select, eachrow
import InvertedIndices: Not
# using CUDA

## Time is the (abstract) type for times.
const Time = Union{Int64,Float64}
const RealTime = Float64
const RealTimeGPU = Float32 

## LogLik is the type for log likelihoods, etc.
## We may at some point need more precision....
const LogLik = Float64
const LogLikGPU = Float32 

# some global options that can be changed at run-time 
mutable struct Config
    usethreads::Bool 
end
const CONFIG = Config(false)
setusethreads!(x::Bool) = begin 
    CONFIG.usethreads = x 
end

include("reshape.jl")

export pomp 
include("pomp.jl")

export times, timezero, obs, init_state, states, coef
include("helpers.jl")

export rinit, rinit!
include("rinit.jl")

export rprocess, rprocess!
include("rprocess.jl")

export euler, discrete_time, onestep
include("plugins.jl")

export rmeasure
include("rmeasure.jl")

export logdmeasure, logdmeasure!
include("logdmeasure.jl")

export simulate,simulate_array
include("simulate.jl")

export pomp, pfilter
include("pfilter.jl")

export melt
include("melt.jl")

# submodule 
include("examples/Examples.jl")

end # module
