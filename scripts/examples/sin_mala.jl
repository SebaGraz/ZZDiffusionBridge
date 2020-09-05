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
iterpol(t, x1, x2, t1, t2) = x2 - (x2-x1)*(t2-t)/(t2-t1)


"""
Return potential function and its gradient
"""
function get_info(X::SinSDE, S::System, ξ::Vector{Float64}, u::Float64, v::Float64, dt::Float64, N::Int64)
        U = deepcopy(norm(ξ)^2/2) #gaussian part
        ∇U = deepcopy(ξ) #gaussian part
        precision = 1
        dt = dt/precision
        t = dt*0.5
        tt = 0:S.T/(2<<S.L):S.T
        XX = fs_expansion(ξ, u, v, S.L, S.T)
        for i in 1:N
                for i2 in 1:precision
                    Xₜ = interpolate(t, XX[i], XX[i+1], tt[i], tt[i+1])
                    #Xₜ = (XX[i]+ XX[i+1])/2
                    U += 0.5*X.α*(X.α*sin(Xₜ)^2 + cos(Xₜ))*dt
                    # to be otpimized the same way as bouncy particle sampler
                    C = 0.5*(X.α^2*sin(2.0*Xₜ) - X.α*sin(Xₜ))*dt
                    for jj in 1:2^(S.L+1)-1
                            if S.ϕ[jj].lb < t < S.ϕ[jj].ub
                                λ = Λ(t, S.ϕ[jj].i, S.ϕ[jj].j, S.T)
                                ∇U[jj] += λ*C
                            end
                    end
                    t += dt
                end
        end
        return ∇U, U
end

function mala_sampler(X::AbstractModel, T::Float64, L::Int64,
                u::Float64, v::Float64, niter::Int64,
                ξ_init =fill(0.0, 2<<L - 1), dt = T/(2<<L))

        ξ₀ = ξ_init
        ξtrace=Array{Float64}(undef, length(ξ_init), niter)
        ξtrace[:,1] = ξ₀
        N = Int(T/dt)
        S = System(L, T, ξ₀, zero(ξ₀))
        @show dt
        @show N
        ∇U₀, U₀ = get_info(X, S, ξ₀, u, v, dt, N)
        τ = 0.015
        count = 0
        for i in 2:niter
                ξ₁  = ξ₀ .- τ*∇U₀ .+ sqrt(2*τ)*randn(length(ξ₀))
                # @assert ξ₁ == ξ₀
                # @show dt
                # @show N
                ∇U₁, U₁ = get_info(X, S, ξ₁, u, v, dt, N)
                # @show i
                # @show U₁
                # @show U₀
                # @assert U₁ == U₀
                # @assert ∇U₁ == ∇U₀
                # @show (norm(ξ₁ - ξ₀ + τ*∇U₀)^2 - norm(ξ₀ - ξ₁ + τ*∇U₁)^2)/(4τ)
                acc_rej =  U₀ - U₁ + (norm(ξ₁ - ξ₀ + τ*∇U₀)^2 - norm(ξ₀ - ξ₁ + τ*∇U₁)^2)/(4τ)
                # @show acc_rej
                # @assert acc_rej == 0.0
                if acc_rej > log(rand())
                        ∇U₀, U₀, ξ₀ = ∇U₁, U₁, ξ₁
                        count += 1
                end
                if i % 100 == 0
                    if count/100 <= 0.7
                        τ = max(0.01, τ - 0.001)
                    else
                        τ = min(1.0, τ + 0.001)
                    end
                    #println("Adaptive step: ar: $(count/100), new tau: $τ")
                    count = 0
                end
                ξtrace[:,i] = ξ₀
        end
        ξtrace
end
function runall_mala(SHORT = false, hist = false)
    Random.seed!(0)
    T = 50.0
    niter = 135000
    L = 6
    α = 0.7#sin
    u = -π
    v = 3π
    X = SinSDE(α, L, T)
    time  = @elapsed (XX = mala_sampler(X, T, L, u, v, niter))
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
    #plot_mala(XX, L, T, u, v)
    time, XX
end

function plot_mala(XX, L, T, u, v)
    tt = 0:T/(2<<L):T
    p = plot(leg = false)
    for i in 1:500:size(XX)[2]
        Y = fs_expansion(XX[:,i], u, v, L, T)
        plot!(tt, Y, alpha = 0.5)
    end
    display(p)
    p
end
#time_mala, XX_mala = runall_mala()
