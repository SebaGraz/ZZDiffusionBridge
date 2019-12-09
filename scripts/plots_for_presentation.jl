include("../src/ZZDiffusionBridge.jl")

using Random



function BM(ξ, L)
    u = 0.0
    v = randn()
    T = 1.0
    x =fs_expansion(ξ, u, v, L, T)
    t = range(0.0, stop=1, length=2<<(L) + 1)
    a = plot(leg=false)
    plot!(a, t, x, label = "random realization of X⁵", alpha = 0.5, color = :red)
    xaxis!(a, "t")
    scatter!(a, t, x, color= :black, markersize=2)
    return a
end

L =5
plot(BM(randn(2<<L-1), L))


function FSbasis(L)
    a = plot()
    plot!(a, 0:0.001:1, 0:0.001:1, label = "first $L levels basis functions", alpha = 0.5)
    for i in 0:L
        for k in 0:2^(i)-1
            plot!(a, 0:0.001:1, [Λ(t, i , k, 1.0) for t in 0:0.001:1], alpha = 0.5, leg=false)
            xaxis!(a, "t")
        end
    end
    return a
end


L = 6
Random.seed!(4)
ξ = randn(2^(L+1) - 1)
plot( FSbasis(L), BM(ξ, L), layout = (1,2))
savefig("../../output/bm_6levels.pdf")

using LinearAlgebra
using Plots
using Random

function zomming()
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
        push!(a, plot(tt[id], xx[id], color = :darkred, legend = false, xaxis = "t"))
        #, label = "zoom X$(2^rep)", legend =:bottomleft ))
        i = Int(ceil((i + f)*0.5))
        dt =  Int(ceil(dt/2))
    end
    return a
end

Random.seed!(5)
a = runall()


plot(plot(a[1], a[3], a[5], a[7], size = (750,500),
       layout = (2,2)))
savefig("../../output/bmzooming.pdf")


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
function vel(x)
    θ = fill(0, length(x[1]))
    θ[1] = 1
    for i in 2:length(x[1])
        θ[i] = -θ[i-1]
    end
    return θ
end

Random.seed!(0)
x = ZigZagSampler(Normal(0.0,1.0), 0.0, 40.0)
θ = vel(x)
a1 = plot(x[1], x[2], leg= false,xaxis = "t", yaxis = "\\xi")
hline!(a1, [0.0], alpha = 0.4, color = :violet)
a2 =  plot(x[1], θ, linetype=:steppost, leg= false, xaxis = "t", yaxis = "\\theta")

plot( a1, a2, layout = (1,2), size = (500, 200))
savefig("../../output/zigzag1.pdf")
