
#Faber Schauder Functions structure
include("faber.jl")

#structures for models and ZigZag sampling schemes
include("types.jl")

#functions for Faber Schauder expansion by finite elements representation
include("fs_expansion.jl")


#ZigZag sampler
include("zz_sampler.jl")

#plotting zigzag sampler
include("plotting.jl")


include("examples/example_sin.jl")
#Experiment
