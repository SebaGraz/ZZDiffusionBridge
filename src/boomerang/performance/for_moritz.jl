using BenchmarkTools
include("../../ZZDiffusionBridge.jl")
include("../../boomerang/exp_1.jl")

"""
Takes a quasi-oredered vector `v` with ordering `p` (ordered from the second element of v).
Update the ordering such that also the first element is ordered as above. This function replaces
findmin() which should be computationally more expensives.
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

"""
For Global Boomerang
Buonds the Poisson rate with a constant rate
"""
function λbar2(ξ::Vector{Float64}, θ::Vector{Float64}, ∇Ubar::Vector{Float64})
    res = 0.0
    for i in 1:length(ξ)
        res += sqrt(ξ[i]*ξ[i] + θ[i]*θ[i])* ∇Ubar[i]
    end
    return res
end

"""
For local Boomerang
Buonds the Poisson rate with a constant rate
"""
function λbar_ind(ξ::Float64, θ::Float64, ∇Ubar_i::Float64 )
    θ_ubs = sqrt(ξ*ξ + θ*θ)
    return θ_ubs*∇Ubar_i
end



"""
    ∇U_tilde_ind!(∇U::Vector{Float64}, ϕ::Vector{Fs}, α::Float64, L::Int64, T::Float64, u::Float64, v::Float64)
for local Boomerang
estimates the ith component of the gradient of the potential `∇U`
"""
function ∇U_tilde_ind!(∇U::Vector{Float64}, i0::Int64, ξ::Vector{Float64}, ϕ::Vector{Fs}, α::Float64, L::Int64, T::Float64, u::Float64, v::Float64)
        t = MCintegration(ϕ[i0])
        XX = fs_expansion(t, ξ, ϕ, u, v, L, T)
        ϕ_t = Λ(ϕ[i0], t)
        ∇U[i0] = 0.5*ϕ[i0].range*ϕ_t*(α*α*sin(2.0*(XX)) - α*sin(XX))
end

"""
    ∇U_tilde!(∇U_tilde::Vector{Float64}, ξ::Vector{Float64}, ϕ::Vector{Fs}, α::Float64, L::Int64, T::Float64, u::Float64, v::Float64)
for Global boomerang
estimates the full gradient of the potential `∇U`
"""
function ∇U_tilde!(∇U_tilde::Vector{Float64}, ξ::Vector{Float64}, ϕ::Vector{Fs}, α::Float64, L::Int64, T::Float64, u::Float64, v::Float64)
    for n in 1:length(ϕ)
        t = MCintegration(ϕ[n])
        XX = fs_expansion(t, ξ, ϕ, u, v, L, T)
        ϕ_t = Λ(ϕ[n], t)
        ∇U_tilde[n] = 0.5*ϕ[n].range*ϕ_t*(α*α*sin(2.0*(XX)) - α*sin(XX))
    end
end

"""
    R(x, v , ∇U_tilde)
for the Global boomeranf
contour reflections
"""
function R(ξ::Vector{Float64}, θ::Vector{Float64} , ∇U_tilde::Vector{Float64})
    return θ - 2*dot(∇U_tilde, θ)/dot(∇U_tilde, ∇U_tilde)*∇U_tilde
end

##### for the local boomerang
##### the reflection on the velocity
##### is like the Zig-Zag:  R(θ_i) = - θ_i
##### Do not need the complicated evaluation above



function local_boomerang(α::Float64, c::Float64, T::Float64, L::Int64, u::Float64, v::Float64, clock::Float64)
    Random.seed!(0)
    N = 2^(L+1)-1
    ξ = zeros(N)
    θ = randn(N)
    t = 0.0
    ϕ = generate(L, T)
    ∇Utilde = zeros(N)
    [∇U_tilde_ind!(∇Utilde, i, ξ, ϕ, α, L, T, u, v) for i in 1:N]
    ∇Ubar = ∇U_bar(α, ϕ, L, T)
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
            ξ, θ =  boomerang_traj(ξ, θ, τ_ref0)
            t += τ_ref0
            θ[i_ref0] = randn()
            push!(Ξ, (Skeleton2(copy(ξ), copy(θ), t)))
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
            ∇U_tilde_ind!(∇Utilde, i0, ξ, ϕ, α, L, T, u, v)
            acc_ratio = max(∇Utilde[i0]*θ[i0], 0)/λ_bar[i0] #max not necessary
            if acc_ratio > rand()
                θ[i0] = -θ[i0]
                push!(Ξ, (Skeleton2(copy(ξ), copy(θ), t)))
                λ_bar[i0] = λbar_ind(ξ[i0], θ[i0], ∇Ubar[i0])
                τ[i0] = event_λ_const(λ_bar[i0])
            else
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



function global_boomerang(α::Float64, c::Float64, T::Float64, L::Int64, u::Float64, v::Float64, clock::Float64)
    Random.seed!(0)
    ξ = zeros(2^(L+1)-1)
    θ = randn(2^(L+1)-1)
    t = 0.0
    ϕ = generate(L, T)
    ∇Utilde = zeros(2^(L+1)-1)
    ∇Ubar = ∇U_bar(α, ϕ, L, T) # does not depend on x, v so precompile
    λ_bar = λbar2(ξ, θ, ∇Ubar)
    τ0 = event_λ_const(λ_bar)
    τ_ref = event_λ_const(c)
    Ξ = [Skeleton2(copy(ξ), copy(θ), t)]
    while t < clock
        if τ_ref < τ0
            ξ, θ =  boomerang_traj(ξ, θ, τ_ref)
            t += τ_ref
            θ = randn(length(θ))
            push!(Ξ, (Skeleton2(copy(ξ), copy(θ), t)))
            λ_bar = λbar2(ξ, θ, ∇Ubar)
            τ0 = event_λ_const(λ_bar)
            τ_ref = event_λ_const(c)
        else
            ξ, θ =  boomerang_traj(ξ, θ, τ0)
            t +=  τ0
            τ_ref -= τ0
            ∇U_tilde!(∇Utilde, ξ, ϕ, α, L, T, u, v)
            acc_ratio = max(dot(∇Utilde, θ), 0)/λ_bar
            if acc_ratio > rand()
                θ[:] = R(ξ, θ, ∇Utilde)
                push!(Ξ, (Skeleton2(copy(ξ), copy(θ), t)))
                λ_bar = λbar2(ξ, θ, ∇Ubar)
                τ0 = event_λ_const(λ_bar)
            else
                τ0 = event_λ_const(λ_bar)
            end
        end
    end
    return Ξ
end


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
