
module Examples 

using ..PartiallyObservedMarkovProcesses

import DataFrames: DataFrame
import CSV: File
import Distributions: Normal, LogNormal, logpdf, Binomial, NegativeBinomial

export brownian_motion
export gompertz
export sir
export rmca
export parus_data

include("parus_data.jl")
include("gompertz.jl")
include("brownian_motion.jl")
include("sir.jl")
include("rmca.jl")

end 
