include("../src/ZZDiffusionBridge.jl")

using Random



function BM(ξ, L)
    u = 0.0
    v = randn()
    T = 1.0
    x =fs_expansion(ξ, u, v, L, T)
    t = range(0.0, stop=1, length=2<<(L) + 1)
    a = plot(leg=false)
    plot!(a, t, x, label = "random realization of X⁵" )
    scatter!(a, t, x, color= "blue", markersize=2)
    return a
end

L =5
plot(BM(randn(2<<L-1), L))


function FSbasis(L)
    a = plot()
    plot!(a, 0:0.001:1, 0:0.001:1, label = "first $L levels basis functions")
    for i in 0:L
        for k in 0:2^(i)-1
            plot!(a, 0:0.001:1, [Λ(t, i , k, 1.0) for t in 0:0.001:1], leg=false)
        end
    end
    return a
end

L = 5
ξ = randn(2^(L+1) - 1)
plot(FSbasis(L), BM(ξ, L), layout = (1,2))
png("./output/bm_5levels")

using LinearAlgebra
using Plots
using Random

function zomming()
    p = plot()
    dt = 10^(-4)
    t = 0.0:dt:100.0
    d = length(t)
    ν = 2.0
    μ = 0.0
    σ = 1.0
    N = randn(d-1)
    x = zeros(d)
    for i = 1:d-1
        x[i+1] = x[i] + ν*(μ -x[i])*dt  + σ*sqrt(dt)*N[i]
    end
    return (t, x)
end

function runall()
    tt, xx = zomming()
    f = length(tt)
    i = 1
    dt = 1000
    a = []
    for rep in 0:6
        id = i:dt:f
        push!(a, plot(tt[id], xx[id], label = "zoom X$(2^rep)", legend =:bottomleft ))
        i = Int(ceil((i + f)*0.5))
        dt =  Int(ceil(dt/2))
        println(i)
        println(f)
        println(dt)
    end
    return a
end

a = runall()
plot(a[1], a[3], a[5], a[7], layout = (2,2))
png("./output/bmzooming.png")


struct Normal
    μ::Float64
    σ::Float64
    v::Float64
    function Normal(μ, σ)
        new(μ, σ, inv(σ))
    end
end


function ZigZagSampler(N::Normal, ξ0::Float64, T::Float64)
    t = 0.0 ; θ = 1; ξ = ξ0
    Ξ = [ξ]
    Δ = [t]
    while (t<T)
        τ = sqrt.(-2*N.σ*log(rand()) + max(0, θ*(ξ - N.μ))^2) - θ*(ξ - N.μ)
        t = t + τ
        ξ = ξ + θ*τ
        push!(Ξ, ξ)
        push!(Δ, t)
        θ = -θ
    end
    return(Δ, Ξ)
end


x = ZigZagSampler(Normal(0.0,1.0), 0.0, 20.0)




function vel(x)
    θ = fill(0, length(x[1]))
    θ[1] = 1
    for i in 2:length(x[1])
        θ[i] = -θ[i-1]
    end
    return θ
end

θ = vel(x)


a1 = plot(x[1], x[2], leg= false)
a2 =  plot(x[1], θ, linetype=:steppost, leg= false)
plot(a1, a2, layout = (1,2), size = (750,250))
png("../output/zigzag1.png")
