module Examples

using ..PartiallyObservedMarkovProcesses

export parus_data
include("parus_data.jl")

export gompertz
include("gompertz.jl")

export brownian_motion
include("brownian_motion.jl")

export sir
include("sir.jl")

export rmca
include("rmca.jl")

export drmca
include("drmca.jl")

end
