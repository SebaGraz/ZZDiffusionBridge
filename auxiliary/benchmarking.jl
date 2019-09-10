using BenchmarkTools
using Random
include("./../scripts/examples/sin.jl")

Random.seed!(1)
function benchmark_sin()
    T = 100.0
    clock = 5000.0
    L = 6
    α = 0.7 #sin
    u = - π
    v = + 3π
    X = SinSDE(α)
    zz_sampler(X, T, L, u, v, clock)
    return
end

@benchmark benchmark_sin()

include("./../scripts/examples/ou.jl")

Random.seed!(1)
function benchmark_ou()
    T = 10.0
    clock = 500.0
    L = 6
    ν = 1.0
    μ = -5.0
    u = -1.0
    v = 2.0
    X = OUSDE(μ, ν, L, T)
    zz_sampler(X, T, L, u, v, clock)
    return
end

@benchmark benchmark_ou()


include("./../scripts/examples/exponential_growth.jl")
function benchmark_expgrowth()
    T = 200.0
    clock = 400.0
    L = 6
    K = 2000
    r = 0.1
    β = 0.1
    u = -log(50)/β     # end points in the lamperti transform
    v= -log(1000)/β
    X = LogGrowthSDE(r, K, β, L)   #end points in the lamperti tranform
    XX = zz_sampler(X, T, L, u, v, clock)
end

@benchmark benchmark_expgrowth()


include(./zzbb.jl)
function benchmark_bb()
    T = 1.0
    clock = 500.0
    L = 6
    u = 0.0
    v = 0.0
    XX = zz_sampler(BB(), T, L, u, v, clock)
end


@benchmark benchmark_bb()
