
#########################################
# Mala for comparison with adaptation   #
#########################################
using Random
include("../../src/faberschauder.jl")
include("../../src/asvar.jl")
include("../../src/pdmp.jl")

function get_info(ϕ, ∇ϕ, ξ, L, α, T)
    U = 0.5*dot(ξ, ξ) #gaussian part
    ∇U = deepcopy(ξ) #gaussian part
    dt = T/(2<<L)
    t = dt*0.5
    tt = 0:T/(2<<L):T
    N = Int(T/dt)
    for i in 1:N
        Xₜ = dotψ(ξ, t, L, T)
        U += ϕ(Xₜ, ξ, t, L, α, T)*dt
        for i in 2:(length(ξ)-1)
            l = lvl(i, L)
            k = (i - 1) ÷ (2 << l)
            δ = T / (1 << (L - l))
            if δ*k < t < δ*(k+1)
                ∇U[i] += ∇ϕ(Xₜ, l, t, ξ,  L, α, T)*dt
            end
        end
        t += dt
    end
    ∇U[1] = ∇U[end] = 0.0
    return ∇U, U
end

function mala_sampler(ϕ, ∇ϕ, ξ₀::Vector{Float64}, niter::Int64, args...)
        ξtrace = Array{Float64}(undef, length(ξ₀), niter)
        ξtrace[:,1] = ξ₀
        ∇U₀, U₀ = get_info(ϕ, ∇ϕ, ξ₀, args...)
        ξ₁ = deepcopy(ξ₀)
        τ = 0.015
        count = 0
        for i in 2:niter
                # initial and final point fixed
                for j in 2:length(ξ₀)-1
                    ξ₁[j]  = ξ₀[j] - τ*∇U₀[j] + sqrt(2*τ)*randn()
                end
                ∇U₁, U₁ = get_info(ϕ, ∇ϕ, ξ₁, args...)
                acc_rej =  U₀ - U₁ + (norm(ξ₁ - ξ₀ + τ*∇U₀)^2 - norm(ξ₀ - ξ₁ + τ*∇U₁)^2)/(4τ)
                # @show acc_rej
                # @assert acc_rej == 0.0
                if acc_rej > log(rand())
                        ∇U₀, U₀, ξ₀ = deepcopy(∇U₁), deepcopy(U₁), deepcopy(ξ₁)
                        count += 1
                end
                if i % 100 == 0
                    if count/100 <= 0.6
                        τ = max(0.0001, τ - 0.0001)
                    else
                        τ = min(1.0, τ + 0.0001)
                    end
                    #println("Adaptive step: ar: $(count/100), new tau: $τ")
                    count = 0
                end
                ξtrace[:,i] = ξ₀
        end
        ξtrace
end

# sin application
function ϕ(Xₜ, ξ, t, L, α, T)
    0.5*α*(α*sin(Xₜ)^2 + cos(Xₜ))
end

function ∇ϕ(Xₜ, l, t, ξ,  L, α, T)
    Λ(t, L - l, T)*0.5*(α^2*sin(2.0*Xₜ) - α*sin(Xₜ))
end

function run_mala(df, niter, α)
    L = 6
    T = 50.0
    n = (2 << L) + 1
    ξ₀ =  randn(n)
    u, v = 0.0, 0.0  # initial and fianl point
    ξ₀[1] = u / sqrt(T)
    ξ₀[end] = v / sqrt(T)
    dt = 1/(2 << L)
    mala_time = @elapsed (out = mala_sampler(ϕ, ∇ϕ, ξ₀, niter, L, α, T))
    burn1 = 1000
    xT_2 = sqrt(T)/2*(out[Int((length(ξ₀)+1)/2), burn1+1:niter])
    nn = length(xT_2)
    ave = cumsum(xT_2)./(Float64.(1:nn))
    for batches in [50]
        ess = ESS(out[:,burn1:end], n_batches = batches)
        push!(df, Dict(:sampler => "MALA", :alpha => α, :T => niter, :nbatches => batches,
            :stat => "ess_mean", :y => sum(ess[2:end-1])/(length(ess)-2), :runtime => mala_time))
        push!(df, Dict(:sampler => "MALA", :alpha => α, :T => niter, :nbatches => batches,
            :stat => "ess_median", :y => median(ess[2:end-1]), :runtime => mala_time))
        push!(df, Dict(:sampler => "MALA", :alpha => α, :T => niter, :nbatches => batches,
            :stat => "ess_min", :y => minimum(ess[2:end-1]), :runtime => mala_time))
        push!(df, Dict(:sampler => "MALA", :alpha => α, :T => niter, :nbatches => batches,
            :stat => "ess_x_T2", :y => ess[Int((length(ξ₀)+1)/2)], :runtime => mala_time))
    end
    # S = T*(0:n)/(n+1)
    # p1 = lines(S, [dotψ(xx[:,end], s, L, T) for s in S], linewidth=0.3)
    # for i in 1:100:size(xx)[2]
    #     lines!(p1, S, [dotψ(xx[:, i], s, L, T) for s in S], linewidth=0.1, alpha = 0.1)
    # end
    # display(p1)
    return df, ave
end

function data_collection_mala(df)
    Random.seed!(0)
    niter = 250000
    mala = []
    for α in [0.1, 0.3, 0.5]
        df, x = run_mala(df, niter, α)
        push!(mala, x)
    end
    return df, mala
end
df, mala = data_collection_mala(df)
# With Automatic Integration (TOO SLOW)
using CSV, JLD2
CSV.write("./scripts/zz_diff_bridges/compare/benchamrk_complete.csv", df)
@save  "./scripts/zz_diff_bridges/compare/covenrgence_mala.jld2" mala
error("Stop here")

cumsum(randn(5,2))


using QuadGK

# sin application
function ϕ(ξ, t, L, α, T)
    Xₜ = dotψ(ξ, t, L, T)
    0.5*α*(α*sin(Xₜ)^2 + cos(Xₜ))
end
function ∇ϕ(l, t, ξ,  L, α, T)
    Xₜ = dotψ(ξ, t, L, T)
    #println(Λ(t, L - l, T)*0.5*(α^2*sin(2.0*Xₜ) - α*sin(Xₜ)))
    0.5*(α^2*sin(2.0*Xₜ) - α*sin(Xₜ))*Λ(t, L - l, T)
end

function get_info(ϕ, ∇ϕ, ξ, L, α, T)
    U = 0.5*dot(ξ, ξ) #gaussian part
    ∇U = deepcopy(ξ) #gaussian part
    U += quadgk(t ->  ϕ(ξ, t, L, α, T), 0.0, T, rtol = 1e-1)[1]
    for i in 2:(length(ξ)-1)
        l = lvl(i, L)
        k = (i - 1) ÷ (2 << l)
        δ = T / (1 << (L - l))
        ∇U[i] += quadgk(t -> ∇ϕ(l, t, ξ, L, α, T), δ*k, δ*(k+1), 1e-1)[1]
    end
    ∇U[1] = ∇U[end] = 0.0
    return ∇U, U
end
runall()
