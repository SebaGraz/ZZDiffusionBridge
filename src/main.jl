include("faber.jl")
include("structures.jl")
include("helper.jl")
include("plotting.jl")
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

function runall(SHORT = false)
    α = 0.7
    σ = 1.0
    L = 7
    T = 100.0
    clock = 5000.0
    u = -π
    v = 3π
    N = 1
    @time x = fullzz(L, α, σ, u, v, T, clock, N)
    if SHORT
        return x
    end
    burning = 10    #burning
    f = clock - 1.0; n = 60
    db = (f-burning)/n
    b =  burning:db:f
    plotmixing(x, b, T, L, u, v)
    return x
end



x = runall(false)
png("./output/plot2")
#vline!(25)
#error("STOP HERE")
