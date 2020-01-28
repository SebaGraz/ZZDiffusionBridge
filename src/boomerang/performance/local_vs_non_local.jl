using Plots
include("../../ZZDiffusionBridge.jl")
include("../../boomerang/exp_3.jl")
include("../../boomerang/exp_2.jl")
include("../../boomerang/exp_1.jl")
"""
Takes a oredered vector `v` with ordering `p` (ordered from the second element
from the smallest value to the biggest).
Update the ordering such that also the first element is ordered as above. This function replaces
findmin() which should be more expensives.
#TODO generaliz for the ith compoment
"""
function first_event_ordering(v::Vector{Float64}, p::Vector{Int64})
    p1 = p[1]
    i_new = searchsortedfirst(v[p[2:end]], v[p1])
    if i_new == 1
        return
    else
        p[1:i_new - 1] = p[2: i_new]
        p[i_new] = p1
    end
end


#
# function event_ordering(ith::Int64, v::Vector{Float64}, p::Vector{Int64})
#     p1 = v[ith]
#     if v[ith]
#     i_new = searchsortedfirst(v[p[1:ith]], v[p1])
#     if i_new == 1
#         return
#     else
#         p[1:i_new - 1] = p[2: i_new]
#         p[i_new] = p1
#     end
# end





function local_boomerang(α::Float64, c::Float64, T::Float64, L::Int64, u::Float64, v::Float64, clock::Float64)
    Random.seed!(0)
    N = 2^(L+1)-1
    ξ = zeros(N)
    θ = randn(N)
    t = 0.0
    ϕ = generate(L, T)
    #initialize quantities
    ∇Utilde = [∇U_tilde_ind(i, ξ, ϕ, α, L, T, u, v) for i in 1:N]
    ∇Ubar = ∇U_bar(α, ϕ, L, T)# does not depend on x, v so precompile
    λ_bar =  [λbar_ind(ξ[i], θ[i], ∇Ubar[i]) for i in 1:N] #vector
    τ = event_λ_const.(λ_bar) #vector
    τ_ref = [event_λ_const(c) for i in 1:N] #vector
    p = sortperm(τ)
    p_ref   = sortperm(τ_ref)
    Ξ = [Skeleton2(copy(ξ), copy(θ), t)]
    τ0, i0 = τ[p[1]], p[1]
    τ_ref0, i_ref0 = τ_ref[p_ref[1]], p_ref[1]
    while t < clock
        if τ_ref0 < τ0
            #println("STEP: refreshment with τ_ref = ", τ_ref)
            ξ, θ =  boomerang_traj(ξ, θ, τ_ref0)
            t += τ_ref0
            θ[i_ref0] = randn()
            push!(Ξ, (Skeleton2(copy(ξ), copy(θ), t)))
            ### RESCALE
            #τ = rescale!(i_ref0, N, τ_ref0, τ_ref)
            for i in 1:(i_ref0 - 1)
                τ_ref[i] -= τ_ref0
                τ[i] -= τ_ref0
            end
            for i in  (i_ref0 + 1): N
                τ_ref[i] -=  τ_ref0
                τ[i] -= τ_ref0
            end
            ## Draw new time
            λ_bar[i_ref0] =  λbar_ind(ξ[i_ref0], θ[i_ref0], ∇Ubar[i_ref0])
            τ[i_ref0] = event_λ_const(λ_bar[i_ref0])
            τ_ref[i_ref0] = event_λ_const(c)
            first_event_ordering(τ_ref, p_ref)
            p = sortperm(τ) #TODO # I have to substitute this step
        else
            ξ, θ =  boomerang_traj(ξ, θ, τ0)
            t +=  τ0
            τ_ref .-= τ0
            ∇Utilde[i0] = ∇U_tilde_ind(i0, ξ, ϕ, α, L, T, u, v)
            acc_ratio = max(∇Utilde[i0]*θ[i0], 0)/λ_bar[i0] #max not necessary
            #if   !(0 <= acc_ratio < 1) #DEBUG
                #println("λ_bar: ", λ_bar)
                #println("lambda tilde ",  max(dot(∇Utilde, θ), 0))
                #println("acc_ratio: ", acc_ratio)
                #error("invalid acc/rej ratio")
            #end
            if acc_ratio > rand()
                #println("STEP: ok accept event time with prob: ", acc_ratio)
                θ[i0] = -θ[i0]
                push!(Ξ, (Skeleton2(copy(ξ), copy(θ), t)))
                λ_bar[i0] = λbar_ind(ξ[i0], θ[i0], ∇Ubar[i0])
                τ[i0] = event_λ_const(λ_bar[i0])
            else
                #println("STEP: ok acc_ratio, reject event time"
                # λ_bar = λbar2(ξ, θ, ∇Ubar) It is not needed
                τ[i0] = event_λ_const(λ_bar[i0])
            end
            for i in 1:(i0-1)
                τ[i] -= τ0
            end
            for i in  (i0 + 1): N
                τ[i] -=  τ0
            end
            first_event_ordering(τ, p)
        end
        τ0, i0 = τ[p[1]], p[1]
        τ_ref0 , i_ref0 =  τ_ref[p_ref[1]], p_ref[1]
    end
    return Ξ
end
#### OK we need just to change  one line ####


function mixed_boomerang(α::Float64, c::Float64, T::Float64, L::Int64, u::Float64, v::Float64, clock::Float64)
    Random.seed!(0)
    N = 2^(L+1)-1
    ξ = zeros(N)
    θ = randn(N)
    t = 0.0
    ϕ = generate(L, T)
    #Q_try = Q(α, ϕ, L, T) #DEBUG
    #if Q_try != Q_try'
    #    error("Q not symmetric")
    #end  # does not depend on x, v so precompile till here
    ∇Utilde = [∇U_tilde_ind(i, ξ, ϕ, α, L, T, u, v) for i in 1:N]
    ∇Ubar = ∇U_bar(α, ϕ, L, T)# does not depend on x, v so precompile
    λ_bar =  [λbar_ind(ξ[i], θ[i], ∇Ubar[i]) for i in 1:N] #vector
    τ = event_λ_const.(λ_bar) #vector
    p = sortperm(τ)
    #println("time: ", τ0) #DEBUG
    #println("a + b*τ0: ", a + b*τ0) #DEBUG
    τ_ref = event_λ_const(c)
    Ξ = [Skeleton2(copy(ξ), copy(θ), t)]
    τ0, i0 = τ[p[1]], p[1]
    while t < clock
        if τ_ref < τ0
            #println("STEP: refreshment with τ_ref = ", τ_ref)
            ξ, θ =  boomerang_traj(ξ, θ, τ_ref)
            t += τ_ref
            θ = randn(N)
            push!(Ξ, (Skeleton2(copy(ξ), copy(θ), t)))
            λ_bar =  [λbar_ind(ξ[i], θ[i], ∇Ubar[i]) for i in 1:N]
            τ = event_λ_const.(λ_bar)
            p = sortperm(τ)
            τ_ref = event_λ_const(c)
        else
            #println("STEP: τ = ", τ0)
            ξ, θ =  boomerang_traj(ξ, θ, τ0)
            t +=  τ0
            τ_ref -= τ0
            ∇Utilde[i0] = ∇U_tilde_ind(i0, ξ, ϕ, α, L, T, u, v) #TODO
            acc_ratio = max(∇Utilde[i0]*θ[i0], 0)/λ_bar[i0] #max not necessary
            #if   !(0 <= acc_ratio < 1) #DEBUG
                #println("λ_bar: ", λ_bar)
                #println("lambda tilde ",  max(dot(∇Utilde, θ), 0))
                #println("acc_ratio: ", acc_ratio)
                #error("invalid acc/rej ratio")
            #end
            if acc_ratio > rand()
                #println("STEP: ok accept event time with prob: ", acc_ratio)
                θ[i0] = -θ[i0]
                push!(Ξ, (Skeleton2(copy(ξ), copy(θ), t)))
                λ_bar[i0] = λbar_ind(ξ[i0], θ[i0], ∇Ubar[i0])
                τ[i0] = event_λ_const(λ_bar[i0])
            else
                #println("STEP: ok acc_ratio, reject event time"
                # λ_bar = λbar2(ξ, θ, ∇Ubar) It is not needed
                τ[i0] = event_λ_const(λ_bar[i0])
            end
            for i in 1:(i0-1)
                τ[i] -= τ0
            end
            for i in  (i0 + 1): N
                τ[i] -=  τ0
            end
            first_event_ordering(τ, p)
        end
        τ0, i0 = τ[p[1]], p[1]
    end
    return Ξ
end


function global_boomerang(α::Float64, c::Float64, T::Float64, L::Int64, u::Float64, v::Float64, clock::Float64)
    Random.seed!(0)
    ξ = zeros(2^(L+1)-1)
    θ = randn(2^(L+1)-1)
    t = 0.0
    ϕ = generate(L, T)
    #Q_try = Q(α, ϕ, L, T) #DEBUG
    #if Q_try != Q_try'
    #    error("Q not symmetric")
    #end  # does not depend on x, v so precompile till here
    ∇Utilde = zeros(2^(L+1)-1)
    ∇Ubar = ∇U_bar(α, ϕ, L, T) # does not depend on x, v so precompile
    λ_bar = λbar2(ξ, θ, ∇Ubar)
    τ0 = event_λ_const(λ_bar)
    #println("time: ", τ0) #DEBUG
    #println("a + b*τ0: ", a + b*τ0) #DEBUG
    τ_ref = event_λ_const(c)
    Ξ = [Skeleton2(copy(ξ), copy(θ), t)]
    while t < clock
        if τ_ref < τ0
            #println("STEP: refreshment")
            ξ[:], θ[:] =  boomerang_traj(ξ, θ, τ_ref)
            t += τ_ref
            θ[:] = randn(length(θ))
            push!(Ξ, (Skeleton2(copy(ξ), copy(θ), t)))
            λ_bar = λbar2(ξ, θ, ∇Ubar)
            τ0 = event_λ_const(λ_bar)
            τ_ref = event_λ_const(c)
        else
            ξ[:], θ[:] =  boomerang_traj(ξ, θ, τ0)
            t +=  τ0
            τ_ref -= τ0
            ∇Utilde[:] = ∇U_tilde(ξ, ϕ, α, L, T, u, v)
            acc_ratio = max(dot(∇Utilde, θ), 0)/λ_bar
            #if   !(0 <= acc_ratio < 1) #DEBUG
                #println("λ_bar: ", λ_bar)
                #println("lambda tilde ",  max(dot(∇Utilde, θ), 0))
                #println("acc_ratio: ", acc_ratio)
                #error("invalid acc/rej ratio")
            #end
            if acc_ratio > rand()
                #println("STEP: ok acc_ratio, accept event time")
                θ[:] = R(ξ, θ, ∇Utilde)
                push!(Ξ, (Skeleton2(copy(ξ), copy(θ), t)))
                λ_bar = λbar2(ξ, θ, ∇Ubar)
                τ0 = event_λ_const(λ_bar)
            else
                #println("STEP: ok acc_ratio, reject event time"
                # λ_bar = λbar2(ξ, θ, ∇Ubar) It is not needed
                τ0 = event_λ_const(λ_bar)
            end
        end
    end
    return Ξ
end



######## MIXING ###########################
function mixing_test()
    α = 0.4
    T = 50.
    L = 6
    c1 = 0.01
    c2 = (2^(L+1)-1)*c1
    u = Float64(-π)
    v = Float64(3π)
    clock  = 20000.
    Ξ1 = local_boomerang(α, c1, T, L, u, v, clock)
    Ξ2 = mixed_boomerang(α, c2, T, L, u, v, clock)
    Ξ3 = global_boomerang(α, c2, T, L, u, v, clock)
    #scaling limits Zig-Zag and bouncy!
end
events_1 = length(Ξ1)
events_2 = length(Ξ2)
events_3 = length(Ξ3)

function find_boomerang_coordinates(Ξ, t)
    tt = [Ξ[i].t for i in 1:length(Ξ)]
    i = searchsortedfirst(tt, t)
    #i = findfirst(x -> x>t ,x in tt)
    ξ, _ = boomerang_traj(Ξ[i-1].ξ, Ξ[i-1].θ, t - Ξ[i-1].t)
    return ξ
end
function plot_boomerang(Ξ, b, L, T, u, v)
    p = Plots.plot(leg = false, colorbar = true)
    dt = range(0.0, stop=T, length=2<<(L) + 1)
    N = length(b)
    P = []
    for i in b
        ξ_interp = find_boomerang_coordinates(Ξ, i)
        dx = fs_expansion(ξ_interp, u, v, L, T)
        push!(P, dx)
    end
    println(length(P) == N)
    p = plot(dt, P, color = :darkrainbow, line_z = (1:N)', linewidth=0.001, alpha = 0.2, leg = false, colorbar = true)
    hline!(p, [n*π for n in -1:2:3], color = :blue)
    display(p)
    return p
end


b = 20:40.0:clock-.1
length(b)
p1 = plot_boomerang(Ξ1, b, L, T, u, v)
p2 = plot_boomerang(Ξ2, b, L, T, u, v)
p3 = plot_boomerang(Ξ3, b, L, T, u, v)

savefig(p1, "./boomerang/performance/output/mix_local.pdf")
savefig(p2, "./boomerang/performance/output/mix_mixed_local.pdf")
savefig(p3, "./boomerang/performance/output/mix_global.pdf")


############################################################

using BenchmarkTools
clock = 1000.
α = 0.4
T = 50.
L = 6
c1 = 0.01
c2 = (2^(L+1)-1)*c1
u = Float64(-π)
v = Float64(3π)


@benchmark local_boomerang($α, $c1, $T, $L, $u, $v, $clock)

# memory estimate:  1.70 GiB
# allocs estimate:  4419733
# --------------
# minimum time:     707.195 ms (12.90% GC)
# median time:      730.503 ms (14.63% GC)
# mean time:        781.007 ms (20.22% GC)
# maximum time:     985.336 ms (35.68% GC)
# --------------
# samples:          7
# evals/sample:     1

@benchmark mixed_boomerang($α, $c2, $T, $L, $u, $v, $clock)

#  memory estimate:  2.02 GiB
#  allocs estimate:  5559079
#  --------------
#  minimum time:     819.092 ms (13.37% GC)
#  median time:      846.517 ms (14.06% GC)
#  mean time:        911.074 ms (20.15% GC)
#  maximum time:     1.114 s (33.83% GC)
#  --------------
#  samples:          6
#  evals/sample:     1

@benchmark global_boomerang($α, $c2, $T, $L, $u, $v, $clock)

#   memory estimate:  34.09 GiB
#   allocs estimate:  96635584
#   --------------
#   minimum time:     28.852 s (8.60% GC)
#   median time:      28.852 s (8.60% GC)
#   mean time:        28.852 s (8.60% GC)
#   maximum time:     28.852 s (8.60% GC)
#   --------------
#   samples:          1
#   evals/sample:     1


#############################################################
