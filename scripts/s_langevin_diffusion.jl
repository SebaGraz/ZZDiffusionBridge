include("../src/ZZDiffusionBridge.jl")
using Distributions
using Plots

struct SinSDE <: AbstractModel
    α::Float64
    # optional precompiled factors
end


"""

`ss` is the stepsize, `n` batchsize for Monte Carlo integration, `nstep' sampler
"""
struct NewSystem
    ξ::Array{Float64, 1}
    ϕ::Vector{Fs}
    L::Int64
    T::Float64
    drift::Array{Float64, 1}
    function NewSystem(L, T)
        new(fill(0.0, 2<<L-1), generate(L, T), L, T, fill(0.0, 2<<L-1))
    end
end


function stoch_langevin(S::NewSystem, X::AbstractModel, ss, nstep, skip, u, v, n)
    d = 2<<S.L - 1
    if ss >= 1
        error("stepsize has to be strictly less the 1")
    end
    Ξ = []
    for i in 1:nstep
        stoch_drift!(S, X, u, v, n)
        S.ξ .=  S.ξ + S.drift*0.5*ss .+ sqrt(ss)*randn(d)  #otherwise: setfield! immutable struct of type NewSystem cannot be changed

        if i%skip == 0
            push!(Ξ, copy(S.ξ))
        end
    end
    Ξ
end





function stoch_drift!(S::NewSystem, X::AbstractModel, u, v, nn)
    for n in 1:2<<S.L-1
        ave = 0.0
        w = 1/nn
        for i in 1:nn
            t = MCintegration(S.ϕ[n])
            XX = fs_expansion(S, t, u, v)
            λ = Λ(t, S.ϕ[n].i , S.ϕ[n].j, S.T)
            ave += w*(0.5*S.ϕ[n].range*λ*(X.α*X.α*sin(2.0*(XX)) - X.α*sin(XX)) + S.ξ[n])
        end
        S.drift[n] = ave
    end
end




function runall(SHORT = false)
    L = 6
    T = 50.0
    S = NewSystem(L, T)
    α = 1.0
    X = SinSDE(α)
    ss = 0.001
    nstep = 10^4
    skip = 10^2
    u = -3.0
    v = 3.0
    n = 5
    Y = stoch_langevin(S, X, ss, nstep, skip, u, v, n)
    if SHORT == false
        burning  = 10
        p = plot(leg = false)
        dt = 0:T/2<<L:T
        for i in burning:length(Y)
            plot!(p, dt, fs_expansion(Y[i], u, v, L, T))
        end
        hline!(p, [n*π for n in -5:2:5])
        display(p)
    end
    return Y, p
end





function fs_expansion(S::NewSystem, t::Float64, u::Float64, v::Float64, n = i -> 2^-(1 + i/2))
        dt = 0:S.T/(2<<S.L):S.T
        k = (searchsortedfirst(dt, t) - 1)
        j0 =  Int(ceil(k/2))-1
        n0 = Faber(S.L, j0)
        if k % 2 != 0
                return interpolate([dt[k], dt[k + 1]], fs_expansion(S.ϕ[n0], S.ξ, u, v, S.L, S.T)[1:2], t)
        else
                return interpolate([dt[k], dt[k + 1]], fs_expansion(S.ϕ[n0], S.ξ, u, v,  S.L, S.T)[2:3], t)
        end
end


Y, p = runall()



error("")
