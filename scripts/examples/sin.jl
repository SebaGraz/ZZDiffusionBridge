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
    δ = S.ϕ[n].δ*(X.α*X.α + X.α)*0.5
    b = S.ξ[n]*S.θ[n]
    if b>0
        return Sol2E(S.θ[n]*S.θ[n]*0.5, b + δ*abs(S.θ[n]), log(ran))
    elseif S.ξ[n]*sign(S.θ[n])*δ <= log(ran) #Case 2: 0 < t < -b/c
        return -log(ran)/(abs(S.θ[n])*δ)
    else    #Case 2: 0 < -b/c < t
        return Sol2E(S.θ[n]*S.θ[n]*0.5, b + δ*abs(S.θ[n]), log(ran)+ S.ξ[n]*S.ξ[n]*0.5)
    end
end



"""
    λratio(n::Int64, S::System, X::SinSDE, u::Float64, v::Float64)

accept reject time drwan from upper bound λbar relative to the coefficient `n`
of model `SinSDE` starting at `u` and ending at `v`
"""
function λratio(n::Int64, S::System, X::SinSDE, u::Float64, v::Float64, t::Float64)
    δ = S.ϕ[n].δ*(X.α*X.α + X.α)*0.5 #always the same, we could save the value
    t = MCintegration(S.ϕ[n])
    XX = fs_expansion(S, t, u, v)
    λ = Λ(t, S.ϕ[n].i, S.ϕ[n].j, S.T)
    acc_rej =  max(0, S.θ[n]*(0.5*S.ϕ[n].range*λ*(X.α*X.α*sin(2.0*(XX)) - X.α*sin(XX)) + S.ξ[n]))/(abs(S.θ[n])*δ + max(0, S.θ[n]*S.ξ[n]))
    # if !(0.0<=acc_rej<=1.0) #DEBUG
    #     println("ratio is :", acc_rej)
    #     error("Ratio outside boundaries")
    # end
    return acc_rej
end

function runall(SHORT = false, hist = false)
    Random.seed!(0)
    T = 50.0
    clock = 10000.0
    L = 6
    α = 0.7#sin
    u = - Float64(π)
    v = 3*Float64(π)
    X = SinSDE(α, L, T)
    XX = zz_sampler(X, T, L, u, v, clock)
    if SHORT == false
        if hist
            burning = 10.0    #burning
            f = clock - 1.0; n = 10000
            db = (f-burning)/n
            b =  burning:db:f
            h = plotmixing(XX, b, T, L, u, v, true)
            #hline!(p, [n*π for n in -3:2:5], color =:blue)
            #display(p)
            return h
        else
            burning = 10.0    #burning
            f = clock - 1.0; n = 200
            db = (f-burning)/n
            b =  burning:db:f
            p = plotmixing(XX, b, T, L, u, v, false)
            hline!(p, [n*π for n in -3:2:5], color =:blue)
            xaxis!(p, "t")
            yaxis!(p, "X_t")
            display(p)
            return XX
        end

    end
end
error("STOP HERE")
x = runall()
savefig("../../output/sin_07.pdf")

#qqplot(Normal, qq[1], label = "level 1", legend = :topleft)
using Distributions
using StatsPlots
using Statistics
histogr(x[:,2])

qqplot(Normal, x[:,6    ])
plot(qqplot(Normal, qq[1], label = "level 0", legend = :topleft, alpha = 0.2), qqplot(Normal, qq[2], label = "level 1", legend = :topleft, alpha = 0.1), qqplot(Normal, qq[3], label = "level 2", legend = :topleft, alpha = 0.2),  size = (500,250),layout = (1,3))
cov(x)
cor(x)
savefig("../../output/qqplots2.pdf")

error("STOP HERE")
using DataFrames
using CSV
df_coef = df = DataFrame(xi_L1=x[:,1] , xi_L2 = x[:,2], xi_L3 = x[:,3],  xi_L4 = x[:,4],  xi_L5 = x[:,5],  xi_L6 = x[:,6])
FILENAME_OUT = "../../output/df_coef.csv"
CSV.write(FILENAME_OUT, df_coef)

cov = DataFrame()
