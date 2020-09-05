
#####################################################
# Try Bouncy particle sampler with exact bounds     #
#####################################################

using ZigZagBoomerang, SparseArrays, LinearAlgebra, Random
#using CairoMakie
include("../../faberschauder.jl")
include("../../../../ICML-boomerang/julia/asvar.jl")
include("../../../../ICML-boomerang/julia/pdmp.jl")

const ZZB = ZigZagBoomerang
using CSV
using DataFrames
using Statistics
# Drift
b(x) = α * sin(x)
# First derivative
b′(x) = α * cos(x)
# Second derivative
b″(x) = -α * sin(x)

struct MyBound
    c::Float64
end
# overloading Poisson times
# WARNING: this function works only for diffusionn bridges and not processes. to generalize it
# change [T^(1.5)/2^((L - lvl(i, L)) * 1.5 + 2) * (α^2 + α) for i in eachindex(x)]
function ZZB.ab(x, θ, c::MyBound, F::BouncyParticle)
    a = c.c + dot([T^(1.5)/2^((L - lvl(i, L)) * 1.5 + 2) * (α^2 + α) for i in eachindex(x)], abs.(θ))
    b1 = dot(x, θ)
    b2 = dot(θ, θ)
    return a, b1, b2
end

function ZZB.adapt!(b::MyBound, x)
    b = MyBound(b.c * x)
end

#Do not reflect the first and the last coefficient
function ZZB.reflect!(∇ϕx, x, θ, F::BouncyParticle)
    θ[2:end-1] .-= (2*dot(∇ϕx[2:end-1], θ[2:end-1])/ZZB.normsq(∇ϕx[2:end-1]))*∇ϕx[2:end-1]
    θ
end
#Do not refresh the first and the last coefficient
function ZZB.refresh!(θ, F::BouncyParticle)
    for i in eachindex(θ)
        θ[i] = randn()
    end
    θ[1] = θ[end] = 0.0
    θ
end
ZZB.sλ(∇ϕx, θ, F) = ZZB.λ(∇ϕx, θ, F)
ZZB.sλ̄((a, b, c)::NTuple{3}, Δt) = a + ZZB.pos(b + c * Δt)


function ∇ϕ!(y, ξ, (L, T); K = 1)
    s = T * (rand())
    x = dotψ(ξ, s, L, T)
    bb = 2b(x) * b′(x) + b″(x)
    for i in eachindex(ξ)
        if i == (2 << L) + 1    # final point
            y[i] = 0.0
        elseif i == 1   # initial point
            y[i] = 0.0
        else
            l = lvl(i, L)
            k = (i - 1) ÷ (2 << l)
            δ = T / (1 << (L - l))
            if δ*k <= s <=   δ*(k + 1)
                y[i] = 0.5 * δ * Λ(s, L - l, T) * bb + ξ[i]
            else
                y[i] =  ξ[i]
            end
        end
    end
    y
end

function run_bps(df, T′)
    c_bps = MyBound(0.0)
    n = (2 << L) + 1
    ξ0, θ0 = randn(n), randn(n)
    θ0[end] = θ0[1] = 0.0 # fix final point
    u, v = 0.0, 0.0  # initial and fianl point
    ξ0[1] = u / sqrt(T)
    ξ0[end] = v / sqrt(T)
    λref_bps = 1.0
    Γ = sparse(1.0I, n, n)
    B = BouncyParticle(Γ, ξ0*0, λref_bps)
    bps_time = @elapsed((bps_trace, uT, acc) = pdmp(∇ϕ!, 0.0, ξ0, θ0, T′, c_bps, B, (L, T); adapt=false, factor=0.0))
    t_bps = getindex.(bps_trace.events, 1)
    x_bps = getindex.(bps_trace.events, 2)
    v_bps = getindex.(bps_trace.events, 3)
    burn1 = 100.0
    i1 = findfirst(x -> x > burn1, t_bps)
    xT_2 = sqrt(T)/2*getindex.(x_bps[i1:end], Int((length(ξ0)+1)/2))
    ave_xT_2 = cumsum(((xT_2[1:end-1] .+ xT_2[2:end])/2).*diff(t_bps[i1:end]))./(t_bps[i1+1:end] .- burn1)
    for batches in [50]
        ess = ess_pdmp_components(t_bps[i1:end] .- burn1, x_bps[i1:end], v_bps[i1:end], n_batches = batches)
        push!(df, Dict(:sampler => "BPS", :alpha => α, :T => T′, :nbatches => batches,
            :stat => "ess_mean", :y => sum(ess[2:end-1])/(length(ess)-2), :runtime => bps_time))
        push!(df, Dict(:sampler => "BPS", :alpha => α, :T => T′, :nbatches => batches,
            :stat => "ess_median", :y => median(ess[2:end-1]), :runtime => bps_time))
        push!(df, Dict(:sampler => "BPS", :alpha => α, :T => T′, :nbatches => batches,
            :stat => "ess_min", :y => minimum(ess[2:end-1]), :runtime => bps_time))
        push!(df, Dict(:sampler => "BPS", :alpha => α, :T => T′, :nbatches => batches,
            :stat => "ess_x_T2", :y => ess[Int((length(ξ0)+1)/2)], :runtime => bps_time))
    end
    return df, (ave_xT_2, t_bps)
end

function data_collection_bps(df)
    Random.seed!(0)
    ave = []
    for α′ in [0.1, 0.3, 0.5]
        global α = α′
        global L = 6
        global T = 50.0
        T′ = 25000
        df, x1 = run_bps(df, T′)
        push!(ave, x1)
    end
    return df, ave
end
df, bps = data_collection_bps(df)
CSV.write("./scripts/zz_diff_bridges/compare/benchamrk_bps_final.csv", df)
@save "convergence_bps.jld2" bps
error("")
# dt = 10.0
# xx = ZZB.trajectory(discretize(out1, dt))
#
# using Makie
# using CairoMakie
# S = T*(0:n)/(n+1)
# p1 = lines(S, [dotψ(xx.x[end], s, L, T) for s in S], linewidth=0.3)
# for ξ in xx.x[1:5:end]
#     lines!(p1, S, [dotψ(ξ, s, L, T) for s in S], linewidth=0.1, alpha = 0.1)
# end
# display(p1)
