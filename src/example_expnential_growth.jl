###################################
####    LOGISTIC GROWTH SDE  ######
###################################
"""
    LogGrowthSDE <: AbstractModel

dX_t = r X_t (1- X_t/K) dt + β X_t dB_t
r := exponential growth
K := saturation parameter
β := multiplicative noise factor
"""
struct LogGrowthSDE <: AbstractModel
    r::Float64
    K::Float64
    β::Float64
    a1::Float64 #precomputing some useful factors
    a2::Float64 #precomputing some useful factors
    function LogGrowthSDE(r, K, β)
        new(r, K, β, 2*r*r/(β*K) , 2*r*r/(β*K*K))
    end
end


# sampling scheme for each stochastic differential equation
sampling_scheme(::LogGrowthSDE) = SubSampling()
# dependence structure for each stochastic differential equation
dependence_strucute(::LogGrowthSDE) = PartialIndependence()


"""
    waitgaus(ξ, θ, u)

Poisson rate  λ(s) = max(0, θ(ξ + θs))
invert function Λ(t) = int_0^t λ(s) ds + ln(u)
"""
function waitgaus( ξ::Float64, θ::Float64, u::Float64)
    b = ξ*θ
    if b>0
        return Sol2E(0.5, b, log(u))
    else
        return Sol2E(0.5, b, log(u) + ξ*ξ*0.5)
    end
end
"""
    waitexp(a, b, u)

invert function  Λ(t) = a/b (e^(bt) - 1) + ln(u)
"""
function waitexp(a, b, u)
    if (log(u)*b/a > 1)
        return error("exponential gave infinity")
    else
        return log(1-log(u)b/a)/b
    end
end
"""
    λbar(n, S::System, X::LogGrowthSDE , u, v)

Poisson time (upper bound) for the coefficient `n`
of model `LogGrowthSDE` starting at `u` and ending at `v`
"""
function λbar(n, S::System, X::LogGrowthSDE , u::Float64, v::Float64, t::Float64)
    w1 = waitgaus(S.ξ[n], S.θ[n], rand())
    S.b1[n] = minimum(fs_expansion(S.ϕ[n], S.ξ, u, v, S.L, S.T))
    S.b2[n] = minimum(fs_expansion(S.ϕ[n], S.θ, u, v, S.L, S.T))
    S.tt[n] = t
    if S.θ[n] > 0
        w2 = waitexp(0.5*S.ϕ[n].δ*X.a1*exp(-X.β*S.b1[n]), -X.β*S.b2[n], rand())
        return min(w1,w2)
    else
        w2 = waitexp(0.5*S.ϕ[n].δ*X.a2*exp(-2X.β*S.b1[n]), -2X.β*S.b2[n], rand())
        return min(w1, w2)
    end
end


"""
    λratio(n::Int64, S::System, X::LogGrowthSDE, u::Float64, v::Float64)

accept reject time drwan from upper bound λbar relative to the coefficient `n`
of model `LogGrowthSDE` starting at `u` and ending at `v`
"""
function λratio(n::Int64, S::System, X::LogGrowthSDE, u::Float64, v::Float64, t::Float64)
    time = t - S.tt[n]
    ω = MCintegration(S.ϕ[n])
    XX = fs_expansion(S, ω, u, v)
    ϕ = Λ(ω, S.ϕ[n].i, S.ϕ[n].j, S.T)   #need to be changed
    num = max(0, S.θ[n]*(0.5*S.ϕ[n].range*ϕ*(X.a1*exp(-X.β*XX) - X.a2*exp(-2X.β*XX)) + S.ξ[n]))
    den = max(0, S.θ[n]*S.ξ[n]) + max(0, 0.5*S.θ[n]*S.ϕ[n].δ*X.a1*exp(-X.β*S.b1[n])*exp(-X.β*S.b2[n]*time)) + max(0, -0.5*S.θ[n]*S.ϕ[n].δ*X.a2*exp(-2X.β*S.b1[n])*exp(-2X.β*S.b2[n]*time))
    return num/den
end


function runall(Short = false)
    T = 200.0
    clock = 200.0
    L = 7
    K = 2000
    r = 0.1
    β = 0.1
    u = -log(50)/β     # end points in the lamperti transform
    v= -log(1000)/β
    X = LogGrowthSDE(r, K, β)   #end points in the lamperti tranform
    XX = zigzagsampler(X, T, L, u, v, clock)
    if SHORT == false
        burning = 10.0    #burning
        f = clock - 1.0; n = 30
        db = (f-burning)/n
        b =  burning:db:f
        p = plotmixing(XX, b, T, L, u, v, x -> exp(- x*β))
        hline!(p, [K])
        display(p)
    end
    return XX
end

x = runall()
error("STOP HERE")
