#################################################################################
# Comparison of Zig-Zag for diffusion bridges with tailored Poisson rates       #
# and with adapted Poisson rate. Reference: https://arxiv.org/abs/2001.05889.   #
#################################################################################

using ZigZagBoomerang, SparseArrays, LinearAlgebra
#using CairoMakie
include("../../src/faberschauder.jl")
include("../../src/asvar.jl")
include("../../src/pdmp.jl")


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


####################################################################
# Overloading Poisson times in order to have tighter upperbounds   #
####################################################################
struct MyBound
    c::Float64
end
function ZZB.adapt!(b::Vector{MyBound}, i, x)
    b[i] = MyBound(b[i].c * x)
end

"""
    poisson_time(a, b, c, u)
Obtaining waiting time for inhomogeneous Poisson Process
with rate of the form λ(t) = a + (b + c*t)^+, where `c`,`a`> 0 ,`b` ∈ R, `u` uniform random variable
"""
function ZZB.poisson_time((a, b, c)::NTuple{3}, u = rand()) # formula (22)
    if b > 0
        return (-(b + a) + sqrt((b + a)^2 - 2.0 * c * log(u))) / c # positive solution of quadratic equation c*0.5 x^2 + (b + a) x + log(u) = 0
    elseif a * b / c <= log(u)
        return -log(u) / a
    else
        return (-(a + b) + sqrt((a + b)^2 - 2.0 * c * (b * b * 0.5 / c + log(u)))) / c    # positive solution of quadratic equation c*0.5 x^2 + (b + a) x + log(u) + b*b*0.5/c = 0
    end
end

ZZB.sλ̄((a, b, c)::NTuple{3}, Δt) = a + ZZB.pos(b + c * Δt)
"""
    abc(G, i, x, θ, c, Flow)
Returns the constant term `a` and linear term `b` when computing the Poisson times
from the upper upper bounding rates λᵢ(t) = max(a + b*t)^2. The factors `a` and `b`
can be function of the current position `x`, velocity `θ`, tuning parameter `c` and
the Graph `G`
"""
function ZZB.ab(G, i, x, θ, c::Vector{MyBound}, F::ZigZag)
    if i == 1
        a = c[i].c + T^(1.5)*0.5*(α^2 + α) * abs(θ[i])  # initial point
        b1 = θ[i]*(x[i] - x[end])
        b2 = θ[i]*(θ[i] - θ[end])
    elseif i == (2 << L) + 1
        a = c[i].c + T^(1.5)*0.5*(α^2 + α) * abs(θ[i])  # final point
        b1 = θ[i]*(x[i] - x[1])
        b2 = θ[i]*(θ[i] - θ[1])
    else
        l = lvl(i, L)
        a = c[i].c + T^(1.5) / 2^((L - l) * 1.5 + 2) * (α^2 + α) * abs(θ[i]) # formula (22)
        b1 = x[i] * θ[i]
        b2 = θ[i] * θ[i]
    end
    return a, b1, b2
end


function ∇ϕmoving_var1(t, ξ, θ, i, t′, F, L, T) # formula (17)
    if i == (2 << L) + 1    # final point
        s = T * (rand())
        x = dotψmoving(t, ξ, θ, t′, s, F, L, T)
        return -b(ξ[end]*sqrt(T))*sqrt(T) + 0.5 * sqrt(T) * s * (2b(x) * b′(x) + b″(x)) + ξ[i] - ξ[1]
    elseif i == 1   # initial point
        s = T * (rand())
        x = dotψmoving(t, ξ, θ, t′, s, F, L, T)
        return 0.5 * T^(1.5) * (1 - s / T) * (2b(x) * b′(x) + b″(x)) + ξ[1]
    else
        l = lvl(i, L)
        k = (i - 1) ÷ (2 << l)
        δ = T / (1 << (L - l))
        res = 0.0
        ave_fact = 1/(l+1)
        for _ in 1:l+1
            s = δ * (k + rand())
            x = dotψmoving(t, ξ, θ, t′, s, F, L, T)
           res += ave_fact * 0.5 * δ * Λ(s, L - l, T) * (2b(x) * b′(x) + b″(x))
       end
       return res + ξ[i]
    end
end


function ∇ϕmoving_var2(t, ξ, θ, i, t′, F, L, T) # formula (17)
    if i == (2 << L) + 1    # final point
        s = T * (rand())
        x = dotψmoving(t, ξ, θ, t′, s, F, L, T)
        return -b(ξ[end]*sqrt(T))*sqrt(T) + 0.5 * sqrt(T) * s * (2b(x) * b′(x) + b″(x)) + ξ[i] - ξ[1]
    elseif i == 1   # initial point
        s = T * (rand())
        x = dotψmoving(t, ξ, θ, t′, s, F, L, T)
        return 0.5 * T^(1.5) * (1 - s / T) * (2b(x) * b′(x) + b″(x)) + ξ[1]
    else
        l = lvl(i, L)
        Δi = T / (1 << L) # subintervals
        k = (i - 1) ÷ (2 << l)
        δ = T / (1 << (L - l))
        res = 0.0
        ave_fact = 1/(l+1)
        for jj in 1:l+1 # number of points proportional to the length
            s = δ * k + Δi*(jj - 1 + rand()) #
            @assert δ * k < s < δ * (k+1)
            x = dotψmoving(t, ξ, θ, t′, s, F, L, T)
           res += ave_fact * 0.5 * δ * Λ(s, L - l, T) * (2b(x) * b′(x) + b″(x))
       end
       return res + ξ[i]
    end
end


##
function run_zz(df, T′)
    println("zigzag var 1")
    n = (2 << L) + 1
    ξ0 = 0randn(n)
    u, v = 0.0, 0.0  # initial and fianl point
    ξ0[1] = u / sqrt(T)
    ξ0[end] = v / sqrt(T)
    θ0 = rand((-1.0, 1.0), n)
    θ0[end] = θ0[1] = 0.0 # fix final point
    Γ = sparse(1.0I, n, n)
    c = [MyBound(0.0) for i in 1:n]
    #Improved algorithm
    zz_time = @elapsed((trace, (t, ξ, θ), (acc, num), c) = spdmp(∇ϕmoving_var1, 0.0, ξ0, θ0, T′, c, ZigZag(Γ, ξ0 * 0),
        SelfMoving(), L, T, adapt = false))
    t_trace = getindex.(trace.events, 1)
    t_trace = [0.0, t_trace...]
    x_trace, v_trace = resize(trace)
    burn1 = 100.0
    i1 = findfirst(x -> x > burn1, t_trace)
    xT_2 = sqrt(T)/2*getindex.(x_trace[i1:end], Int((length(ξ0)+1)/2))
    ave_xT_2_var = cumsum(((xT_2[1:end-1] + xT_2[2:end])/2).*diff(t_trace[i1:end]))./(t_trace[i1+1:end] .- burn1)
    for batches in [50]
        ess = ess_pdmp_components(t_trace[i1:end] .- burn1, x_trace[i1:end], v_trace[i1:end], n_batches = batches)
        push!(df, Dict(:sampler => "ZZ_var", :alpha => α, :T => T′, :nbatches => batches,
            :stat => "ess_mean", :y => sum(ess[2:end-1])/(length(ess)-2), :runtime => zz_time))
        push!(df, Dict(:sampler => "ZZ_var", :alpha => α, :T => T′, :nbatches => batches,
            :stat => "ess_median", :y => median(ess[2:end-1]), :runtime => zz_time))
        push!(df, Dict(:sampler => "ZZ_var", :alpha => α, :T => T′, :nbatches => batches,
            :stat => "ess_min", :y => minimum(ess[2:end-1]), :runtime => zz_time))
        push!(df, Dict(:sampler => "ZZ_var", :alpha => α, :T => T′, :nbatches => batches,
            :stat => "ess_x_T2", :y => ess[Int((length(ξ0)+1)/2)], :runtime => zz_time))
    end
    println("standard zigzag")
    n = (2 << L) + 1
    ξ0 = 0randn(n)
    u, v = 0.0, 0.0  # initial and fianl point
    ξ0[1] = u / sqrt(T)
    ξ0[end] = v / sqrt(T)
    θ0 = rand((-1.0, 1.0), n)
    θ0[end] = θ0[1] = 0.0 # fix final point
    Γ = sparse(1.0I, n, n)
    c = [MyBound(0.0) for i in 1:n]
    #Improved algorithm
    zz_time = @elapsed((trace, (t, ξ, θ), (acc, num), c) = spdmp(∇ϕmoving, 0.0, ξ0, θ0, T′, c, ZigZag(Γ, ξ0 * 0),
        SelfMoving(), L, T, adapt = false))
    t_trace1 = getindex.(trace.events, 1)
    t_trace1 = [0.0, t_trace1...]
    x_trace, v_trace = resize(trace)
    i1 = findfirst(x -> x > burn1, t_trace1)
    xT_2 = sqrt(T)/2*getindex.(x_trace[i1:end], Int((length(ξ0)+1)/2))
    ave_xT_2 = cumsum(((xT_2[1:end-1] + xT_2[2:end])/2).*diff(t_trace1[i1:end]))./(t_trace1[i1+1:end] .- burn1)
    for batches in [50]
        ess = ess_pdmp_components(t_trace1[i1:end] .- burn1, x_trace[i1:end], v_trace[i1:end], n_batches = batches)
        push!(df, Dict(:sampler => "ZZ", :alpha => α, :T => T′, :nbatches => batches,
            :stat => "ess_mean", :y => sum(ess[2:end-1])/(length(ess)-2), :runtime => zz_time))
        push!(df, Dict(:sampler => "ZZ", :alpha => α, :T => T′, :nbatches => batches,
            :stat => "ess_median", :y => median(ess[2:end-1]), :runtime => zz_time))
        push!(df, Dict(:sampler => "ZZ", :alpha => α, :T => T′, :nbatches => batches,
            :stat => "ess_min", :y => minimum(ess[2:end-1]), :runtime => zz_time))
        push!(df, Dict(:sampler => "ZZ", :alpha => α, :T => T′, :nbatches => batches,
            :stat => "ess_x_T2", :y => ess[Int((length(ξ0)+1)/2)], :runtime => zz_time))
    end
    println("zigzag var 2")
    n = (2 << L) + 1
    ξ0 = 0randn(n)
    u, v = 0.0, 0.0  # initial and fianl point
    ξ0[1] = u / sqrt(T)
    ξ0[end] = v / sqrt(T)
    θ0 = rand((-1.0, 1.0), n)
    θ0[end] = θ0[1] = 0.0 # fix final point
    Γ = sparse(1.0I, n, n)
    c = [MyBound(0.0) for i in 1:n]
    #Improved algorithm
    zz_time = @elapsed((trace, (t, ξ, θ), (acc, num), c) = spdmp(∇ϕmoving_var2, 0.0, ξ0, θ0, T′, c, ZigZag(Γ, ξ0 * 0),
        SelfMoving(), L, T, adapt = false))
    t_trace3 = getindex.(trace.events, 1)
    t_trace3 = [0.0, t_trace3...]
    x_trace3, v_trace3 = resize(trace)
    i1 = findfirst(x -> x > burn1, t_trace3)
    xT_3 = sqrt(T)/2*getindex.(x_trace3[i1:end], Int((length(ξ0)+1)/2))
    ave_xT_3 = cumsum(((xT_3[1:end-1] + xT_3[2:end])/2).*diff(t_trace3[i1:end]))./(t_trace3[i1+1:end] .- burn1)
    for batches in [50]
        ess = ess_pdmp_components(t_trace3[i1:end] .- burn1, x_trace3[i1:end], v_trace3[i1:end], n_batches = batches)
        push!(df, Dict(:sampler => "ZZ_var2", :alpha => α, :T => T′, :nbatches => batches,
            :stat => "ess_mean", :y => sum(ess[2:end-1])/(length(ess)-2), :runtime => zz_time))
        push!(df, Dict(:sampler => "ZZ_var2", :alpha => α, :T => T′, :nbatches => batches,
            :stat => "ess_median", :y => median(ess[2:end-1]), :runtime => zz_time))
        push!(df, Dict(:sampler => "ZZ_var2", :alpha => α, :T => T′, :nbatches => batches,
            :stat => "ess_min", :y => minimum(ess[2:end-1]), :runtime => zz_time))
        push!(df, Dict(:sampler => "ZZ_var2", :alpha => α, :T => T′, :nbatches => batches,
            :stat => "ess_x_T2", :y => ess[Int((length(ξ0)+1)/2)], :runtime => zz_time))
    end
    return df, (ave_xT_2, t_trace1), (ave_xT_2_var, t_trace), (ave_xT_3, t_trace3)
end


function resize(zz_trace)
    x0 = zz_trace.x0
    v0 = zz_trace.θ0
    x_trace = [deepcopy(x0)]
    v_trace = [deepcopy(v0)]
    for i in eachindex(zz_trace.events)
        if i == 1
            x0 = x0 + v0*zz_trace.events[i][1]
        else
            x0 = x0 + v0*(zz_trace.events[i][1] - zz_trace.events[i-1][1])
        end
        v0[zz_trace.events[i][2]] = zz_trace.events[i][4]
        x0[zz_trace.events[i][2]] = zz_trace.events[i][3]
        # if !(v0[zz_trace.events[i][2]] == zz_trace.events[i][4])
        #     println("vel: $(v0[zz_trace.events[i][2]]) vs $(zz_trace.events[i][4])")
        #     error("")
        # end
        # if !(|x0[zz_trace.events[i][2]]| zz_trace.events[i][4])
        #     println("pos: $(x0[zz_trace.events[i][2]]) vs $(zz_trace.events[i][3])")
        #     error("")
        # end
        push!(x_trace, deepcopy(x0))
        push!(v_trace, deepcopy(v0))
    end
    x_trace, v_trace
end

using DataFrames, Random
df = DataFrame(sampler = String[], alpha = Float64[], T = Float64[],
            nbatches = Int64[], stat = String[], y = Float64[], runtime = Float64[])
function data_collection(df)
    Random.seed!(0)
    T′ = 25000.0
    zz_var = []
    zz_var1 = []
    zz = []
    for α′ in [0.1, 0.3, 0.5]
        global α = α′
        global L = 6
        global T = 50.0
        df, x1, x2, x3 = run_zz(df, T′)
        push!(zz, x1)
        push!(zz_var, x2)
        push!(zz_var1, x3)
    end
    return df, zz, zz_var, zz_var1
end
using CSV, JLD2
df, zz, zz_var1, zz_var2 = data_collection(df)
#CSV.write("./scripts/zz_diff_bridges/compare/benchamrk_zz_final.csv", df)
@save "./scripts/compare/covenrgence_zz.jld2" zz
@save "./scripts/compare/covenrgence_zz_var.jld2" zz_var1
@save "./scripts/compare/covenrgence_zz_var1.jld2" zz_var2
error("")

#
# using Makie
# using CairoMakie
#
# T′ = 1000
# α =0.3
# n = (2 << L) + 1
# ξ0 = 0randn(n)
# u, v = 0.0, 0.0  # initial and fianl point
# ξ0[1] = u / sqrt(T)
# ξ0[end] = v / sqrt(T)
# θ0 = rand((-1.0, 1.0), n)
# θ0[end] = θ0[1] = 0.0 # fix final point
# Γ = sparse(1.0I, n, n)
# c = [MyBound(0.0) for i in 1:n]
# zz_time = @elapsed((trace, (t, ξ, θ), (acc, num), c) = spdmp(∇ϕmoving, 0.0, ξ0, θ0, T′, c, ZigZag(Γ, ξ0 * 0),
#     SelfMoving(), L, T, adapt = false))
# t_trace = getindex.(trace.events, 1)
# t_trace = [0.0, t_trace...]
# x_trace, v_trace = resize(trace)
# ess = ess_pdmp_components(t_trace, x_trace, v_trace, n_batches = 50)
# #plot
# ts, ξs = splitpairs(discretize(trace, T′/n))
# S = T*(0:n)/(n+1)
# p1 = lines(S, [dotψ(ξ, s, L, T) for s in S], linewidth=0.3)
# for ξ in ξs[1:5:end]
#     lines!(p1, S, [dotψ(ξ, s, L, T) for s in S], linewidth=0.3)
# end
# display(p1)
