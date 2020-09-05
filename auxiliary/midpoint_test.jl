using StatsPlots
using Distributions
using KernelDensity

include("../scripts/examples/ou.jl")
include("../src/plotting.jl")
#### WARNING: u,v,L has to coincide with the one in ../script/examples/ou.jl inside runall
T = 10.0
ν = 1.0
μ = -5.0
L = 6
u = -1.0
v = 2.0

runall()
function empirical_sample(dt, burn, u, v, L, T)
    X = runall(true)
    tt = burn:dt:(X[end].t-0.01)
    Ξ = []
    for i in tt
         push!(Ξ, FindCoordinates(X, i))
    end
    return Ξ
end

Ξ = empirical_sample(0.1, 10.0, u, v, L, T)


function midpoint_sample(Ξ, u, v, T)
    X = []
    for j in 1:length(Ξ)
        y = u/2 + v/2 +  Λ(T/2,0,0,T)* Ξ[j].ξ[1]
        push!(X, y)
    end
    return Array{Float64,1}(X)
end


midp = midpoint_sample(Ξ, u, v, T)
#a = histogram(midp)
x = kde(midp)
aa = plot(x.x,x.density, leg = false, size = (300, 300) )


#### Real density
# dX_t = ν(μ - X_t) dt + dB_t
#step
a = -5.0
b = -1.0
T2 = T/2
function β(t , a, b, v, T)
    a + 2b*exp(b*(T-t))*(v + a/b*(1 - exp(b*(T-t))))/(exp(2b*(T-t))-1)
end

function int_B(t ,a, b, T2, step, T)
    som = 0
    for i in t:step:(T2 - step)
        som += step*(b - 2b*exp(2b*(T - (i + step/2)))/(exp(2b* (T - (i + step/2))) - 1))
    end
    return som
end

function μ_(a, b, T2, step, u, v, T)
    int1 = 0.0
    for t in 0:step:T2
        int1 += step*exp(int_B(t + step/2, a, b, T2, step, T))*β(t + step/2, a, b, v, T)
    end
    int1 += u*exp(int_B(0, a, b,  T2, step, T))
end

function σ_(T2, step1, step2, a, b, T)
    int1 = 0.0
    for t in 0:step1:T2
        int1 += step1*exp(2*int_B(t + step1/2, a, b, T2, step2, T))
    end
    return int1
end

#mean
μ = μ_(a, b, T2, 0.001, u, v, T)

#variance
σ = σ_(T2, 0.001, 0.001, a, b, T)





plot!(aa, Normal(μ, sqrt(σ)), label = "true distribution" )
savefig("./output/midpoint300.pdf")
