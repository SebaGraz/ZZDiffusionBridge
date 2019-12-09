include("../../src/ZZDiffusionBridge.jl")
using LinearAlgebra
"""
    OUSDE <: AbstractModel

dX_t = ν(μ - X_t)dt + dB_t
ν := intensity
μ := mean reversion
M::Array{Float64, 2} := matrix whose element i,j is equal to ∫_0^T ϕ_i ϕ_j dt
V::Vector{Float64} :=  vector whose element i is equal to ∫_0^T ϕ_i dt
bound1::Vector{Float64} :=  vector whose element i is equal to ∫_0^T upbar{ϕ}_i ϕ_j dt
bound2::Vector{Float64} :=  vector whose element i is equal to ∫_0^T downbar{ϕ}_i ϕ_j dt
"""
struct OUSDE <: AbstractModel
    μ::Float64
    ν::Float64
    M::Array{Float64, 2}
    V::Vector{Float64}
    bound1::Vector{Float64}
    bound2::Vector{Float64}
    function OUSDE(μ, ν, L, T)
        new(μ, ν, generate_matrix(L, T), generate_vector(L, T), generate_bound1(L,T), generate_bound2(L,T))
    end
end

dependence_strucute(::OUSDE) = PartialIndependence()
sampling_scheme(::OUSDE) = Regular()



#Poisson rates
"""
    wait_gengaus(a,b,u)
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

"""
    λbar(n, S::System, X::OUSDE , u, v)

Poisson time (real bound) for the coefficient `n`
of model `OUSDE` starting at `u` and ending at `v`

λbar = λ No upper bound and no accept reject step
"""
function λbar(n, S::System, X::OUSDE , u::Float64, v::Float64, t::Float64)
    a = S.θ[n]*(X.ν*X.ν*(dot(X.M[n,:], S.ξ) + X.bound1[n]*v + X.bound2[n]*u - X.μ*X.V[n]) + S.ξ[n])
    b = S.θ[n]*((dot(X.M[:,n], S.θ))*X.ν*X.ν + S.θ[n])
    return wait_gengaus(a, b, rand())
end

####
function runall(SHORT = false)
    Random.seed!(1)
    T = 10.0
    clock = 1000.0
    L = 6
    ν = 1.0
    μ = -5
    u = -1.0
    v = 2.0
    X = OUSDE(μ, ν, L, T)
    XX = zz_sampler(X, T, L, u, v, clock)
    if SHORT == false
        burning = 10.0    #burning
        f = clock - 1.0; n = 100
        db = (f-burning)/n
        b =  burning:db:f
        p = plotmixing(XX, b, T, L, u, v)
        hline!(p, [-5.0], color = :blue)
        xaxis!(p, "t")
        yaxis!(p, "X_t")
        display(p)
        #plot the mean of the process
    end
    return XX
end

error("STOP HERE")

runall()
savefig("../../output/ou_1_m5.pdf")
