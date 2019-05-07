include("faber.jl")
include("structures.jl")
include("helper.jl")
include("plotting.jl")
using Random



#n = 2^i - 1
L = 8


function BM(ξ, L)
    a = plot(0:0.001:1,[cies(t, ξ, L, 1.0, 0.0, 0.0) for t in 0:0.001:1], leg= false)
    tt = range(0,1,length= 2^(L+1) + 1)
    scatter!(a, tt, [cies(t, ξ, L, 1.0, 0.0, 0.0) for t in tt], color= "blue", markersize=2)
    return a
end


plot(BM(ξ, 8))


function FSbasis(L)
    a = plot()
    for i in 0:L
        for k in 0:2^(i)-1
            plot!(a, 0:0.001:1, [Λ(t, i , k, 1.0) for t in 0:0.001:1], leg=false)
        end
    end
    return a
end


#plot(FSbasis(1))
i = L
ξ = randn(2^(L+1) - 1)

i = 8
plot(FSbasis(i), BM(ξ, i), layout = (1,2))
png("./output/bm8")

using LinearAlgebra
using Plots
using Random

p = plot()
dt = 0.001
t = 0:dt:10
d = length(t)
b = 1.0
σ = 1.0
for j in 1:1
    N = randn(d-1)
    global x = zeros(d)
    for i = 1:d-1
        x[i+1] = x[i] + b*(sin(x[i]))*dt  + σ*sqrt(dt)*N[i]
    end
    plot!(p, t, x, color= RGB(j/10, 0.5, 1-j/10), linewidth=0.3, alpha = 1.0)
end
plot!(p, leg= false)
display(p)
x[1:1001]
plot(0:dt:1, x[1:1001])


x[1:101]
plot(0:dt:0.1, x[1:101])



d = length(t)-1
Δt = t[end] - t[end-1]
p = plot()
for j in 1:10
    dW = randn(d)
    W = cumsum(dW)
    Wend = W[end]
    W = vcat(0.0, W)*sqrt(Δt)
    plot!(p, W)
end
plot!(p, leg=false)
display(p)
