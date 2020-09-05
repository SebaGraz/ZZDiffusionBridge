include("../../src/ZZDiffusionBridge.jl")
using LinearAlgebra
"""
    SinSDE <: AbstractModel

dX_t = α sin(X_t) + dB_t
α := attraction intensity
"""
struct SinSDE <: AbstractModel
    α::Float64
    V::Vector{Float64}
    # optional precompiled factors
    function SinSDE(α, L, T)
        new(α, generate_vector(L, T))
    end
end

"""
waiting times from poisson rate of the form λ(t) = a + (b + ct)^+
with c,a > 0, b ∈ ℝ
"""
function waiting_time(a,b,c)
    ran = rand()
    if b>0
        return Sol2E(c*0.5, b + a, log(ran))
    elseif a*b/c <= log(ran)
        return -log(ran)/a
    else
        Sol2E(c*0.5, a + b, log(ran) + b*b*0.5/c)
    end
end

"""
        λbar(S::System, X::SinSDE, u, v)
Bouncy particle sampler optimized for sparse operations
Poisson time (upper bound) of model `SinSDE`
starting at `u` and ending at `v`

invert function Λ(t) = at + int_0^t max(0, b+cs) ds + ln(ran)
with a = dot(abs.(θ),δ),   b = dot(ξ,θ), c = dot(θ,θ)
"""
function λbar(S::System, X::SinSDE, u::Float64, v::Float64)
    δ = [S.ϕ[n].δ*(X.α*X.α + X.α)*0.5 for n in 1:length(S.θ)]
    a = dot(abs.(S.θ), δ) #>0
    b = dot(S.ξ, S.θ)
    c = dot(S.θ, S.θ)
    waiting_time(a,b,c)
end

function stochastic_gradient(S::System, X::SinSDE, u::Float64, v::Float64, n_mc::Int64)
    Utilde = deepcopy(S.ξ)
    for i in 1:n_mc
        t = MCintegration(S.ϕ[1])
        XX = fs_expansion(S, t, u, v)
        C = (X.α*X.α*sin(2.0*(XX)) - X.α*sin(XX))
        for jj in 1:2^(S.L+1) - 1
            if S.ϕ[jj].lb < t < S.ϕ[jj].ub
                λ = Λ(t, S.ϕ[jj].i, S.ϕ[jj].j, S.T)
                Utilde[jj] += (0.5*S.ϕ[jj].range*λ*C)/n_mc
            end
        end
    end
    Utilde
end

"""
    λratio(n::Int64, S::System, X::SinSDE, u::Float64, v::Float64)
BOUCY PARTICLE SAMPLER
accept reject time drwan from upper bound λbar relative to the coefficient `n`
of model `SinSDE` starting at `u` and ending at `v`
"""
function λratio(S::System, X::SinSDE, Utilde::Vector{Float64})
    δ = [S.ϕ[i].δ*(X.α*X.α + X.α)*0.5 for i in 1:length(S.θ)] #always the same, we could save the value
    a = dot(abs.(S.θ),δ)
    b = dot(S.ξ, S.θ)
    acc_rej =  max(0, dot(S.θ, Utilde))/(a + max(0, b))
      if !(0.0<=acc_rej<=1.0) #DEBUG
          println("ratio is :", acc_rej)
          error("Ratio outside boundaries")
      end
    return acc_rej
end

"""
Bouncy reflection
"""
function reflect!(S::System, Utilde::Vector{Float64})
    θ_new = S.θ - 2*(dot(Utilde,S.θ)/dot(Utilde,Utilde))*Utilde
    S.θ .= θ_new
end


λref(λref::Float64) = -log(rand())/λref



"""
Bp_sampler
"""
function bp_sampler(X::AbstractModel, T::Float64, L::Int64, u::Float64, v::Float64, clock::Float64, ξ =fill(0.0, 2<<L - 1) , θ = randn(2<<L - 1); λ_ref = 1.0)
    #initialization
    t = 0.0
    S = System(L, T, ξ, θ)
    τ0 = 0.0
    n0 = 0
    Ξ = Skeleton[]
    τ0 = λbar(S, X, u, v)
    τref = λref(λ_ref)
    n_mc = 1
    while t < clock
        if τ0 < τref
            S.ξ .+=  S.θ*τ0
            t += τ0
            τref -= τ0
            Utilde = stochastic_gradient(S, X, u, v, n_mc)
            if λratio(S, X, Utilde) > rand()
                reflect!(S, Utilde)
                push!(Ξ, (Skeleton(copy(S.ξ), t)))
                τ0 = λbar(S, X, u, v)
            else
                τ0 = λbar(S, X, u, v)
            end
        else
            S.ξ .+=  S.θ*τref
            t += τref
            randn!(S.θ)
            τ0 = λbar(S, X, u, v)
            τref = λref(λ_ref)
        end
    end
    return Ξ
end


function runall_bouncy(SHORT = false, hist = false)
    Random.seed!(0)
    T = 50.0
    #clock = 150.0
    clock = 8000.0
    L = 6
    α = 0.7#sin
    u = - Float64(π)
    v = 3*Float64(π)
    X = SinSDE(α, L, T)
    time = @elapsed (XX = bp_sampler(X, T, L, u, v, clock))
    # if SHORT == false
    #     if hist
    #         burning = 10.0    #burning
    #         f = clock - 1.0; n = 10000
    #         db = (f-burning)/n
    #         b =  burning:db:f
    #         h = plotmixing(XX, b, T, L, u, v, true)
    #         #hline!(p, [n*π for n in -3:2:5], color =:blue)
    #         #display(p)
    #         return h
    #     else
    #         burning = 10.0    #burning
    #         f = clock - 1.0; n = 200
    #         db = (f-burning)/n
    #         b =  burning:db:f
    #         p = plotmixing(XX, b, T, L, u, v, false)
    #         hline!(p, [n*π for n in -3:2:5], color =:blue)
    #         xaxis!(p, "t")
    #         yaxis!(p, "X_t")
    #         display(p)
    #         return XX
    #     end
    # end
    time, XX
end
#time_bouncy, XX_bouncy = runall_bouncy()
#time_bouncy, XX_bouncy = runall_bouncy()
