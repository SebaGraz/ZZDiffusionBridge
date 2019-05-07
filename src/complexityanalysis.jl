include("faber.jl")
include("structures.jl")
include("helper.jl")
include("plotting.jl")
using BenchmarkTools
using StatsBase

function fullzz(L, α, σ, u, v, T, clock , N=1)
    d =  2^(L+1) - 1 ; t = 0.0 ;
    k = [GFSbase(i, L, T) for i in 1:d] #inizialize basis
    ξ = zeros(d)
    θ = fill(1, d)
    Ξ = Skeleton[]
    τ = zeros(d)
    while t < clock
        τ0, i0 = findmin(τ)
        ξ = ξ + θ*τ0
        t = t + τ0
        if λratio(k[i0], ξ, θ, τ0, α, σ, u, v, N) > rand()
            push!(Ξ, Skeleton(t,ξ))
            θ[i0] = -θ[i0]
        end
        τ[i0] = TimeAbs(k[i0], ξ[i0], θ[i0], α, σ, rand())
        for j in 1: i0 - 1          #can be improved without if statement
            if j in k[i0].nhb
                τ[j] = TimeAbs(k[j], ξ[j], θ[j], α, σ, rand())
            else
                τ[j] = τ[j] - τ0
            end
        end
        for j in i0 + 1:d           #can be improved without if statement
            if j in k[i0].nhb
                τ[j] = TimeAbs(k[j], ξ[j], θ[j], α, σ, rand())
            else
                τ[j] = τ[j] - τ0
            end
        end
    end
    return(Ξ)
end

function runallcomplexity(L)
    α = 0.7
    σ = 1.0
    T = 100.0
    clock = 100.0
    u = -π
    v = 3π
    N = 1
    println(" Time for Level = ", L)
    @benchmark fullzz(L, α, σ, u, v, T, clock, N)

end


runallcomplexity(1)

for L in (100:100:1100)
    runall(L)
end

L=8
α = 0.7
σ = 1.0
T = 100.0
clock = 100.0
u = -π
v = 3π
N = 1
@benchmark fullzz(L, α, σ, u, v, T, clock, N)


#L = 1 33.852 ms (9.19% GC)
#L = 2 75.048 ms (7.19% GC)
#L = 4 305.282 ms (6.74% GC)
#L = 6 1.411 s (4.30% GC)
#L = 8 8.400 s (3.66% GC
