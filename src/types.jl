#TYPES
"""
    AbstractModel

abstract type for models
"""
abstract type AbstractModel end

"""
    SinSDE <: AbstractModel

dX_t = α sin(X_t) + dB_t
α := attraction intensity
"""
struct SinSDE <: AbstractModel
    α::Float64
    # precompiled factors
end

"""
    LogGrowthSDE <: AbstractModel

dX_t = r X_t (1- X_t/K) dt + β X_t dB_t
r := exponential growth
K := saturation parameter
β := multiplicative noise factor
"""
struct LogGrowthSDE <: AbstractModel
    r::Float64
    K::Float64
    β::Float64
    a1::Float64 #precomputing some useful factors
    a2::Float64 #precomputing some useful factors
    function LogGrowthSDE(r, K, β)
        new(r, K, β, 2*r*r/(β*K) , 2*r*r/(β*K*K))
    end
end

"""
    OUSDE <: AbstractModel

dX_t = ν(μ - X_t)dt + dB_t
ν := intensity
μ := mean reversion
"""
struct OUSDE <: AbstractModel
    μ::Float64
    ν::Float64
end
"""
    AbstractDependenceStructure

Types inheriting from abstract type `AbstractDependenceStructure`
"""
abstract type AbstractDependenceStructure end

"""
    FullIndependence <: AbstractDependenceStructure

Type acting as a Flag for Full independent Zig-Zag sampler
"""
struct FullIndependence <: AbstractDependenceStructure end

"""

Type acting as a Flag for partial independent Zig-Zag sampler
"""
struct PartialIndependence <: AbstractDependenceStructure end

"""
    SamplingScheme

Abstact type for Sampling scheme
"""
abstract type SamplingScheme end

"""
    SubSampling <: SamplingScheme

If you cannot sample from the inhomogeneous Poisson rate, sample it by
subsampling from Poisson rate with higher intensity (bound) and accept reject the event.
"""
struct SubSampling <: SamplingScheme end

"""
    Regular <: SamplingScheme

For linear sdes, where the actual imhogeneous Poisson rate can be sampled directly.
"""
struct Regular <: SamplingScheme end

"""
    System

contains all the information needed for the ZigZag sampler
    ξ::Vector{Float64} := vector for the position of the coefficients
    θ::Vector{Float64} := vector containing the velocities (need to be changes in float64)
    ϕ::Vector{Fs} := Faber Schauder functions information (see Fs)
    τ::Vector{Float64} := vector containing waiting time (could be sorted)
    L::Int64 := Number of Levels
    T::Float64 := time length of the bridge
    b1::Vector{Float64} := free vector needed for linear growth sde
    b2::Vector{Float64} := free vector needed for linear growth sde
    tt::Vector{Float64} := free vector needed for linear growth sde
"""
struct System
    ξ::Vector{Float64} #
    θ::Vector{Float64}
    ϕ::Vector{Fs}
    τ::Vector{Float64}
    L::Int64
    T::Float64
    b1::Vector{Float64} #bad programming
    b2::Vector{Float64} #bad programming
    tt::Vector{Float64} #bad programming
    function System(L::Int64, T::Float64, ξ = fill(0.0, 2<<L - 1), θ = fill(1.0, 2<<L - 1), b1 = fill(0.0, 2<<L - 1), b2 = fill(0.0, 2<<L - 1))
        new(ξ, θ, generate(L, T), fill(0.0, 2<<L - 1), L, T, b1, b2, fill(0.0, 2<<L - 1))
    end
end


"""
    Skeleton
container for output zigzag sampler
"""
struct Skeleton
    ξ::Vector{Float64}
    t::Float64
end
