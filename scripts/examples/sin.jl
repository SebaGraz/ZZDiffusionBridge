include("../../src/ZZDiffusionBridge.jl")

"""
    SinSDE <: AbstractModel

dX_t = α sin(X_t) + dB_t
α := attraction intensity
"""
struct SinSDE <: AbstractModel
    α::Float64
    # optional precompiled factors
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
    δ = S.ϕ[n].δ*(X.α*X.α + X.α)*0.5   #always the same, we could save the value
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
    δ = S.ϕ[n].δ*(X.α*X.α + X.α)*0.5 #always the same, we could save the value
    t = MCintegration(S.ϕ[n])
    XX = fs_expansion(S, t, u, v)
    λ = Λ(t, S.ϕ[n].i , S.ϕ[n].j, S.T)
    return max(0, S.θ[n]*(0.5*S.ϕ[n].range*λ*(X.α*X.α*sin(2.0*(XX)) - X.α*sin(XX)) + S.ξ[n]))/(abs(S.θ[n])*δ + max(0, S.θ[n]*S.ξ[n]))
end


function runall(SHORT = false)
    T = 100.0
    clock = 10000.0
    L = 6
    α = 0.7 #sin
    u = - π
    v = + 3π
    X = SinSDE(α)
    XX = zz_sampler(X, T, L, u, v, clock)
    if SHORT == false
        burning = 10.0    #burning
        f = clock - 1.0; n = 50
        db = (f-burning)/n
        b =  burning:db:f
        p = plotmixing(XX, b, T, L, u, v)
        hline!(p, [n*π for n in -5:2:5])
        display(p)
    end
    return XX
end

x = runall()

error("STOP HERE")
png("output/sin_07.png")
