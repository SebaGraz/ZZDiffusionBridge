using BenchmarkTools
include("../src/ZZDiffusionBridge.jl")
struct BB <: AbstractModel
    # optional precompiled factors
end

dependence_strucute(::BB) = FullIndependence()
sampling_scheme(::BB) = Regular()

#Poisson rates
"""
obtaining waiting time for Inhomogeneous Poisson Process
with rate of the form λ(t) = (a + b*t)^+, `a`,`b` ∈ R, `u` random variable
"""
function wait_gengaus(a,b,u)
    if b > 0
        if a < 0
            τ = sqrt(-log(u)*2.0/b) - a/b
        else #a[i]>0
            τ = sqrt((a/b)^2 - log(u)*2.0/b) - a/b
        end
    elseif  b == 0
        if a > 0
            τ = -log(u)/a
        else #a[i] <= 0
            τ = Inf
        end
    else #b[i] < 0
        if a <= 0
            τ = Inf
        elseif -log(u) <= -a^2/b + a^2/(2*b)
            τ = - sqrt((a/b)^2 - log(u)*2.0/b) - a/b
        else
            τ = Inf
        end
    end
end


function λbar(n, S::System, ::BB , u::Float64, v::Float64, t::Float64)
    return wait_gengaus(S.ξ[n]*S.θ[n], 1.0, rand())
end

#This function will not be necessary after fixing issue #4
function update_events!(n::Int64, τ0::Float64, S::System, X::AbstractModel, ::Regular, ::FullIndependence, acc::Bool, u::Float64, v::Float64, t::Float64)
    S.τ[n] = λbar(n, S, X , u, v, t)
    for i in 1:(n-1)
        S.τ[i] -= τ0
    end
    for i in  (n + 1): length(S.ϕ)
        S.τ[i] -=  τ0
    end
end

using Random
Random.seed!(1)
function runall(SHORT = false)
    T = 1.0
    clock = 500.0
    L = 8
    u = 0.0
    v = 0.0
    @time XX = zz_sampler(BB(), T, L, u, v, clock)
    if SHORT == false
        burning = 10.0    #burning
        f = clock - 1.0; n = 200
        db = (f-burning)/n
        b =  burning:db:f
        p = plotmixing(XX, b, T, L, u, v)
        display(p)
        #plot the mean of the process
    end
end




error("STOP HERE")

runall(false)

png("../output/bm00.png")
