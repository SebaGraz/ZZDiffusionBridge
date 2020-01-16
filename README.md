# ZZDiffusionBridge
ZigZag sampler is used to explore the conditional measure of a diffusion process. This measure lays in a high dimensional space (infinite dimensional) and could differ significantly from a Gaussian measure (which in our case is the reference measure). The repo is currently under development. To see and edit the notes see https://www.overleaf.com/2848134699gmmwzrpsyvsx .


## Overview
The main function introduced (available [here](src/zz_sampler.jl)) is  

```julia
zz_sampler(X::AbstractModel, T::Float64, L::Int64, u::Float64, v::Float64, clock::Float64; θ = fill(1.0, 2<<L - 1))
```

It takes a diffusion model with its parameters, its initial and final points `u,v`, the truncation level of the infinite summation `L`, the final time of the ZigZag sampler `clock` and as optional input, the velocity vector for each coordinate. The function returns the skeleton of the ZigZag process with the real event times. This is given in the form of `Array{Skeleton, 1}`. 

### Structures
Three main types are defined ([here](src/types.jl)): 
1. the abstract type `::AbstractModel`, inheriting any specific diffusion model
2. the abstract type `::AbstractDependenceStructure` inheriting the types of dependences (for our application we have now two subtypes `::FullIndependence` and `::PartialIndependence`) acting as a flag
3. the abstract type `SamplingScheme` inheriting the types of sampling scheme (`::Regular` and `::PartialIndependence` ) acting as a flag
4. the type `System`, container of all the attributes necessary for the sampler

## Examples
We have three examples right now ([here](/scripts/examples) running the ZigZag sampler for three different SDEs. **how to run them:** easy, just download the repo and run the script in the folder src/examples 

## Developing the ZigZag for a new SDE
If you want to develop your own SDE, you need to:
1. define your  `::ModelName:<AbstractModel` containing the parameters (and optionally some parameters/vectors/matrices needed to the sampler which is better to precompute.)
2. define the function `λbar(n, S::System, X::ModelName , u, v)` and `λratio(n::Int64, S::System, X::ModelName, u::Float64, v::Float64, t::Float64)` (`λratio` only if the sampling scheme is `SubSampling`)
3. assign the flags for the sampler depending wheather the sampler is `Fullindependence` or `PartialIndependence`; `Regular` or `SubSampling`. This is done for example via the flag-funcitons

dependence_strucute(::ModelName) = FullIndependence()
sampling_scheme(::ModelName) = SubSampling()

## Faber-Schauder functions
The file [faber.jl](src/faber.jl) and [fs_expansion.jl](src/fs_expansion.jl) contains all the functions necessary to work with the Faber-Schauder functions and change of basis to finite elements. 

## Tuning the velocities
The script [tune_velocities.jl](scripts/tune_velocities.jl) computes the sample mean of the path integral for the sde dX_t = \alpha sin(X_t)dt + dB_t starting from -3\pi and ending at 3\pi at time 200. By symmetry of the sde, the real mean should be 0 and therefore we can compute the distance between sample and real mean, as a function of the time of the ZigZag sampler. We repeated the experiment for different values of \beta which determines the decay of the velocties 2^(-i\beta) where i is the level. Finally we plot the results. **N.B** this is still unstable and does not bring the desired results. 

## Comparisons
We would like to compare with the **Stochastic gradient Langevin dynamics** [for info here!](https://en.wikipedia.org/wiki/Stochastic_gradient_Langevin_dynamics). This is because, as the methodology we propose, this sampling method use a Monte Carlo estimator for the Gradient of the energy function. The script is [here](scripts/s_langevin_diffusion.jl), you just need to run it. You can play around with the parameters: step size `step_size`, number of steps `nstep`, number of skipped jumps before saving `skip`.  The method is simple to implement and the results do not seem that bad, although you need to tune the discretization step. If it is too little, it does not mix, if it is too large, it diverges from the real density. From here maybe it is worth to implement a Mala (or stochastic mala).

## Results!

Sine SDE with alpha = 0.7

![temp](output/sin_07.png)

Exponential growth model with r = 0.1, b=0.1 K = 2000

![temp](output/exp_growth_01_01_2000.png)

Orstein-Uhlenbeck process with \mu = -5 and \nu= 1.0

![temp](output/ou_1_m5.png)

Stochastic gradient Langevin dynamic for the sine SDE with alpha = 0.7, batch size = 2^(L-1) - 1 

![temp](output/langevin_sampler_sin_07.png)

