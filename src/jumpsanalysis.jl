# COPY AND PASTE OF MAIN SCRIPT
#adding in the main function a vector that keeps track of the number of jumps for each coordinate



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
    jumps = zeros(d)
    while t < clock
        τ0, i0 = findmin(τ)
        jumps[i0] += 1
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
    return(jumps)
end

function runalljumps()
    α = 0.0
    σ = 1.0
    L = 7
    T = 100.0
    clock = 2000.0
    u = -π
    v = 3π
    N = 1
    @time jumps = fullzz(L, α, σ, u, v, T, clock, N)
    return jumps
end



#vline!(25)
#error("STOP HERE")


function jumps(x)
    L=7
    n =zeros(L+1)
    k = 1
    ii = 1
    for i in 0:L
        for j in 0:2^i-1
            n[ii] += x[k]
            k += 1
        end
        ii +=1
    end
    return n
end


x = runalljumps()
jum = jumps(x)
scatter(x)
png("./output/jumpsXcoordinate0.0e1.0")

2^(7+1)
