include("../../src/ZZDiffusionBridge.jl")

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


# dependence structure for each stochastic differential equation
dependence_strucute(::SinSDE) = FullIndependence()
# sampling scheme for each stochastic differential equation
sampling_scheme(::SinSDE) = SubSampling()


"""
        λbar(n, S::System, X::SinSDE, u, v)

Poisson time (upper bound) for the coefficient `n`
of model `SinSDE` starting at `u` and ending at `v`

invert function Λ(t) = at + int_0^t max(0, b+cs) ds + ln(ran)
with a = |θ|*δ,   b = ξ*θ, c = θ^2
"""
function λbar(n, S::System, X::SinSDE, u, v, t::Float64)
    ran = rand()
    #δ = S.ϕ[n].δ*(X.α*X.α + X.α)*0.5
    δ = X.V[n]*(X.α*X.α + X.α)*0.5   #always the same, we could save the value
    b = S.ξ[n]*S.θ[n]
    if b>0
        return Sol2E(S.θ[n]*S.θ[n]*0.5, b + δ*abs(S.θ[n]), log(ran))
    elseif S.ξ[n]*sign(S.θ[n])*δ <= log(ran) #Case 2: 0 < t < -b/c
        return -log(ran)/(abs(S.θ[n])*δ)
    else    #Case 2: 0 < -b/c < t
        return Sol2E(S.θ[n]*S.θ[n]*0.5, b+δ*abs(S.θ[n]), log(ran)+ S.ξ[n]*S.ξ[n]*0.5)
    end
end


"""
    λratio(n::Int64, S::System, X::SinSDE, u::Float64, v::Float64)

accept reject time drwan from upper bound λbar relative to the coefficient `n`
of model `SinSDE` starting at `u` and ending at `v`
"""
function λratio(n::Int64, S::System, X::SinSDE, u::Float64, v::Float64, t::Float64)
    #δ = S.ϕ[n].δ*(X.α*X.α + X.α)*0.5 #always the same, we could save the value
    δ = X.V[n]*(X.α*X.α + X.α)*0.5
    t = MCintegration(S.ϕ[n])
    XX = fs_expansion(S, t, u, v)
    λ = Λ(t, S.ϕ[n].i, S.ϕ[n].j, S.T)
    return max(0, S.θ[n]*(0.5*S.ϕ[n].range*λ*(X.α*X.α*sin(2.0*(XX)) - X.α*sin(XX)) + S.ξ[n]))/(abs(S.θ[n])*δ + max(0, S.θ[n]*S.ξ[n]))
end

function runall(SHORT = false)
    Random.seed!(0)
    T = 50.0
    clock = 10000.0
    L = 6
    α = 0.7#sin
    u = - Float64(π)
    v = + 3*Float64(π)
    X = SinSDE(α, L, T)
    XX = zz_sampler(X, T, L, u, v, clock)
    if SHORT == false
        burning = 10.0    #burning
        f = clock - 1.0; n = 200
        db = (f-burning)/n
        b =  burning:db:f
        p = plotmixing(XX, b, T, L, u, v)
        hline!(p, [n*π for n in -3:2:5], color =:blue)
        xaxis!(p, "t")
        yaxis!(p, "X_t")
        display(p)
    end
    return XX
end
error("STOP HERE")
x = runall()
savefig("../../output/sin_07.pdf")

function qqplot_coef(xx, levels, burn)
    y = Vector{Float64}[]
    T = xx[end].t - 1
    time = burn:1.0:T
    for ii in levels
        n = Faber(ii, 0)
        x = []
        for i in time
            push!(x, FindCoordinates(xx, i).ξ[n])
        end
        push!(y, copy(x))
    end
    return y
end

levels = [0,1,2]
burn = 100
qq = qqplot_coef(x, levels, burn)
#qqplot(Normal, qq[1], label = "level 1", legend = :topleft)
using Distributions
using StatsPlots
plot(qqplot(Normal, qq[1], label = "level 0", legend = :topleft, alpha = 0.2), qqplot(Normal, qq[2], label = "level 1", legend = :topleft, alpha = 0.1), qqplot(Normal, qq[3], label = "level 2", legend = :topleft, alpha = 0.2),  size = (500,250),layout = (1,3))

savefig("../../output/qqplots2.pdf")

error("STOP HERE")


l = @layout [a{.1h};grid(1,3)]
plot(
       plot(annotation=(0.5,0.5, "Q-Qplots"), framestyle = :none),
       qqplot(Normal, qq[1], label = "level 0", legend = :topleft),
       qqplot(Normal, qq[2], label = "level 3", legend = :topleft),
       qqplot(Normal, qq[3], label = "level 5", legend = :topleft),
       layout = l)
