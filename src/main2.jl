#using Ciesielski
using LinearAlgebra
#using Plots
using SparseArrays
include("faber.jl")
include("structures.jl")
include("gaussiankernel.jl")
include("plotting.jl")
#include("Models/Gaussian.jl")
#include("Faber1.jl")
#include("PDSamplers.jl")

#L = levels of the Faber-Schauder basis
#a1,b1 parameters of the sde dX_t = (a + b X_t)dt + dW_t
#T timer of the ZigZagSampler before exiting
#ξ0 Initial points





function ZigZagExperiment( ξ0::Array{Float64, 1}, a1, b1, L, T::Float64)
    B = sparse(FillFingerMatrix(L, a1, b1))
    A = vector1(L, a1, b1)
    d = length(A) ; t = 0.0 ; θ = ones(d);
    ξ = ξ0
    Ξ = [Skeleton(t,ξ)]
    τ = zeros(d)
    while t<T
        a = θ.*(-A - 2(B*ξ))
        b =-2θ.*(B*θ)
        τ = GetTime.(a,b,rand(d))
        τ1, i0 = findmin(τ)
        #k, kk = Faber(i0)
        t = t + τ1
        ξ = ξ + θ*τ1
        push!(Ξ, Skeleton(t, ξ))
        θ[i0] = -θ[i0]
    end
    return Ξ
end



#experiment
function runall(SHORT = false)
    clock = 5000.0
    L = 7
    d = 2^(L+1) - 1 #dimensions given by the level L
    a = 10.0
    b = -1.0
    ZZ = ZigZagExperiment(zeros(d), a, b, L, clock)
    if SHORT
        return ZZ
    end
    #plotting
    burning = 10    #burning
    f = clock - 1.0; n = 60
    db = (f-burning)/n
    b =  burning:db:f
    plotmixing(ZZ, b, 1.0 , L, 0.0, 0.0)
    return ZZ
end






x = runall()
png("./output/OU_10_1.png")


##### Jumps

function jumpsOU(ξ0::Array{Float64, 1}, a1, b1, L, T::Float64)
    B = sparse(FillFingerMatrix(L, a1, b1))
    A = vector1(L, a1, b1)
    d = length(A) ; t = 0.0 ; θ = ones(d);
    ξ = ξ0
    jumps = zeros(d)
    τ = zeros(d)
    while t<T
        a = θ.*(-A - 2(B*ξ))
        b =-2θ.*(B*θ)
        τ = GetTime.(a,b,rand(d))
        τ1, i0 = findmin(τ)
        jumps[i0] += 1
        t = t + τ1
        ξ = ξ + θ*τ1
        θ[i0] = -θ[i0]
    end
    return jumps
end


function runalljumps()
    clock = 5000.0
    L = 7
    d = 2^(L+1) - 1 #dimensions given by the level L
    a = 10.0
    b = -1.0
    jumps = jumpsOU(zeros(d), a, b, L, clock)
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
scatter(x)
png("./output/jumpsXcoordinateOU_10_1.png")
jum = jumps(x)
scatter(jum)
png("./output/jumpsXlevelsOU_10_1.png")
