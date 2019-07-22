# ZZDiffusionBridge
ZigZag sampler is used to explore the conditional measure of a diffusion process. This measure lays an high dimensional space (infinite dimensional) and could differ significantly from the Gaussian measure (which in our case is the reference measure). The repo is currently under developmemt. To see and edit the notes see https://www.overleaf.com/2848134699gmmwzrpsyvsx .


## Overview
the main function introduced (available [here](src/zz_sampler.jl)) is  

```julia
zz_sampler(X::AbstractModel, T::Float64, L::Int64, u::Float64, v::Float64, clock::Float64; θ = fill(1.0, 2<<L - 1))
```

It takes a diffusion Model with its parameters, its initial and final points `u,v`, the truncation level of the infinite summation `L`, the final time of the ZigZag sampler `clock` and as optional input, the velocity vector for each coordinate. The function returns the skeleton of the ZigZag process with the real event times. This is given in the form of `Array{Skeleton, 1}`. 

### Structures
Three main types are defined ([here](src/types.jl)): 
1. the asbtarct type `::AbstractModel`, inheriting a specific diffusion model
2. the abstract type `::AbstractDependenceStructure` inheriting the type of dependences (for our application we have now two subtypes `::FullIndependence` and `::PartialIndependence`) acting a a flag
3. the abstract type `SamplingScheme` inheriting the type of sampling scheme (`::Regular` and `::PartialIndependence` ) acting as a flag
4. the type `System`, container of all the attributes necessary for the sampler
see the issue [#11](issues/11) for changing system. It contains stuff that may be not used for certain models (ionefficient). Maybe this should be a tuple or we should move all the arguments which are model spcific as attributes of `::Model:<AbstractModel`.

## Examples
We have three examples right now. They are [here](/src/examples) and implement the ZigZag sampler for three different sdes.

## Developing the ZigZag for a new SDE
If you want to develop your own sde, you need to:
1. define your  `::ModelName:<AbstractModel` containing the parameters
2. define the function `λbar(n, S::System, X::ModelName , u, v)` and `λratio(n::Int64, S::System, X::ModelName, u::Float64, v::Float64, t::Float64)` (`λratio` only if the sampling scheme is `SubSampling`)
3. assign the flags for the sampler depending wheather the sampler is `Fullindependence` or `PartialIndependence`; `Regular` or `SubSampling`. This is done for example via the flag-funcitons

dependence_strucute(::ModelName) = FullIndependence()
sampling_scheme(::ModelName) = SubSampling()
## Faber-Schauder functions
the file [faber.jl](src/faber.jl) and (fs_expansion.jl) contains all the functions necessary to work with the Faber Schauder functions and change of basis to finite element basis. 

## Tuning the velocities
TODO

