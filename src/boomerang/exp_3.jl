include("../ZZDiffusionBridge.jl")
using LinearAlgebra
using JLD
# in order to recostruct the path we need time and velocity at event time
"""
        struct Skeleton2
            ξ::Vector{Float64} #does saving take time?
            θ::Vector{Float64}
            t::Float64
        end

structure for saving the output of the process. `ξ` position, `θ` velocity, `t` time
"""
struct Skeleton2
    ξ::Vector{Float64} #does saving take time?
    θ::Vector{Float64}
    t::Float64
end
"""
    boomerang_traj(ξ, θ, t)

trajectories which in this case are circles in (ξ, θ)
"""
function boomerang_traj(ξ::Vector{Float64}, θ::Vector{Float64}, t::Float64)
    ξ_new = ξ.*cos(t) + θ.*sin(t)
    θ = θ.*cos(t) - ξ.*sin(t)
    return ξ_new, θ
end


"""
     ∇U_bar(α::Float64, ϕ::Vector{Fs}, L::Int64, T::Float64)

computes upper bound of the gradient of the potential function(-log(density)).
Does not depend on x --> pre-compile
"""
function ∇U_bar(α::Float64, ϕ::Vector{Fs}, L::Int64, T::Float64)
    return [0.5*ϕ[n].δ*(α^2 + α) for n in 1:length(ϕ)]
end


"""
         fs_expansion(t::Float64, ξ::Vector{Float64}, ϕ::Vector{Fs}, u::Float64, v::Float64, L::Int64, T::Float64, n = i -> 2^-(1 + i/2))

Local expansion: find value of the process (Piecewise linear) for any t ∈ [0, T]
does not evaluate the  whole path, but just the points needed for t
"""
function fs_expansion(t::Float64, ξ::Vector{Float64}, ϕ::Vector{Fs}, u::Float64, v::Float64, L::Int64, T::Float64, n = i -> 2^-(1 + i/2))
        dt = 0:T/(2<<L):T
        k = (searchsortedfirst(dt, t) - 1)
        j0 =  Int(ceil(k/2))-1
        n0 = Faber(L, j0)
        if k % 2 != 0
                return interpolate([dt[k], dt[k + 1]], fs_expansion(ϕ[n0], ξ, u, v, L, T)[1:2], t)
        else
                return interpolate([dt[k], dt[k + 1]], fs_expansion(ϕ[n0], ξ, u, v,  L, T)[2:3], t)
        end
end



#Poisson rates
"""
    wait_linear(a,b,u)
obtaining waiting time for Inhomogeneous Poisson Process
with rate of the form λ(t) = (a + b*t)^+, `a`,`b` ∈ R, `u` random variable
"""
function wait_linear(a::Float64, b::Float64, u::Float64)
    if b > 0
        if a < 0
            τ = sqrt(-log(u)*2.0/b) - a/b
        else #a[i]>0
            τ = sqrt((a/b)^2 - log(u)*2.0/b) - a/b
        end
    elseif  b == 0
        if a > 0
            τ = -log(u)/a
        else #a[i] <= 0
            τ = Inf
        end
    else #b[i] < 0
        if a <= 0
            τ = Inf
        elseif -log(u) <= -a^2/b + a^2/(2*b)
            τ = - sqrt((a/b)^2 - log(u)*2.0/b) - a/b
        else
            τ = Inf
        end
    end
end

"""
    event_λref(c::Float64)
Draw waiting time from homogeneous Poisson rate c. Used for refreshments
"""
function event_λ_const(c::Float64)
    return -log(rand())/c
end

"""
    event_λbar(a::Float64, b::Float64)
Draw waiting time from rate λ(t) = max(0, a + b t)
"""
function event_λbar(a::Float64, b::Float64)
    wait_linear(a, b, rand())
end


# """
#
# Try new bounds, constant
# """
# function λbar2(ξ::Vector{Float64}, θ::Vector{Float64}, ∇Ubar::Vector{Float64})
#     θ_ubs = sqrt.(ξ.*ξ + θ.*θ)
#     return dot(θ_ubs, ∇Ubar)
# end
"""
try new bounds
"""
function λbar_ind(ξ::Float64, θ::Float64, ∇Ubar_i::Float64 )
    θ_ubs = sqrt(ξ*ξ + θ*θ)
    return θ_ubs*∇Ubar_i
end




"""
Try new bound
"""
function acc_rej2(∇U_tilde::Vector{Float64}, global_bound::Float64, θ::Vector{Float64})
    max(0, ∇U_tilde*θ)/global_bound
end

"""
    Λ(ϕ_n::Fs, t::Float64)
fast way to evaluate a basis function `ϕ_n` at a point t (inside its support)
"""
function  Λ(ϕ_n::Fs, t::Float64)
    ϕm = (ϕ_n.lb + ϕ_n.ub)*0.5
    if t < ϕm
        (t - ϕ_n.lb)*ϕ_n.supr/(ϕm - ϕ_n.lb)
    else
        (ϕ_n.ub - t)*ϕ_n.supr/(ϕ_n.ub - ϕm)
    end
end


# """
#     ∇U_tilde(ξ::Vector{Float64}, ϕ::Vector{Fs}, α::Float64, L::Int64, T::Float64, u::Float64, v::Float64)
#
# accept reject time drwan from upper bound λbar relative to the coefficient `n`
# of model `SinSDE` starting at `u` and ending at `v`
# """
# function ∇U_tilde(ξ::Vector{Float64}, ϕ::Vector{Fs}, α::Float64, L::Int64, T::Float64, u::Float64, v::Float64)
#     ∇U_tilde = zeros(length(ϕ))
#     for n in 1:length(ϕ)
#         t = MCintegration(ϕ[n])
#         XX = fs_expansion(t, ξ, ϕ, u, v, L, T)
#         ϕ_t = Λ(ϕ[n], t)
#         ∇U_tilde[n] = 0.5*ϕ[n].range*ϕ_t*(α*α*sin(2.0*(XX)) - α*sin(XX))
#     end
#     return ∇U_tilde
# end

"""
    ∇U_tilde_ind(ξ::Vector{Float64}, ϕ::Vector{Fs}, α::Float64, L::Int64, T::Float64, u::Float64, v::Float64)

accept reject time drwan from upper bound λbar relative to the coefficient `n`
of model `SinSDE` starting at `u` and ending at `v`
"""
function ∇U_tilde_ind!(∇U::Vector{Float64}, i0::Int64, ξ::Vector{Float64}, ϕ::Vector{Fs}, α::Float64, L::Int64, T::Float64, u::Float64, v::Float64)
        t = MCintegration(ϕ[i0])
        XX = fs_expansion(t, ξ, ϕ, u, v, L, T)
        ϕ_t = Λ(ϕ[i0], t)
        ∇U[i0] = 0.5*ϕ[i0].range*ϕ_t*(α*α*sin(2.0*(XX)) - α*sin(XX))
end



function acc_rej(∇U_tilde::Vector{Float64}, θ::Vector{Float64})
    max(0, ∇U_tilde*θ)/(a + b*τ)
end



# """
#     R_ind(x, v , ∇U_tilde)
#
# contour reflections
# """
# function R_ind(i, θ::Array{Float64}, ∇U_tilde::Array{Float64})
#     return θ[i] - 2*dot(∇U_tilde, θ)/dot(∇U_tilde, ∇U_tilde)*∇U_tilde[i]
# end

"""
    R_ind_2(x, v , ∇U_tilde)

contour reflections
"""
function R_ind_2(i, θ::Float64)
    return -θ
end

function rescale!(i0::Int64, N::Int64, τ0::Float64, τ::Vector{Float64})
    for i in 1:(i0-1)
        τ[i] -= τ0
    end
    for i in  (i0 + 1): N
        τ[i] -=  τ0
    end
    return τ
end

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


function boomerang_ind_refresh(T::Float64, L::Int64, u::Float64, v::Float64, clock::Float64)
    clearconsole()
    N = 2^(L+1)-1
    ξ = zeros(N)
    θ = randn(N)
    t = 0.0
    α = 1.0
    c = 0.05
    ϕ = generate(L, T)
    #initialize quantities
    ∇Utilde = zeros(N)
    [∇U_tilde_ind!(∇Utilde, i, ξ, ϕ, α, L, T, u, v) for i in 1:N]
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
            end
            for i in  (i_ref0 + 1): N
                τ_ref[i] -=  τ_ref0
            end
            τ .-= τ_ref0
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
error("STOP HERE")

Random.seed!(0)
T = 50.0
u = Float64(-π)
v =  Float64(3π)
L = 6
clock = 20000.0
Ξ = boomerang_ind_refresh(T, L, u, v, clock)


save("ouput.jld", "output", Ξ)


output = load("ouput.jld")
Ξ = output["output"]


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
    p = plot(dt, P, color = :darkrainbow, line_z = (1:N)', linewidth=0.001, alpha = 0.2, leg = false, colorbar = true)
    hline!(p, [n*π for n in -1:2:3], color = :blue)
    display(p)
    return p
end



using Plots
b = 19:20.0:clock-.1
length(b)
plot_boomerang(Ξ, b, L, T, u, v)


savefig("boomerang/sin_boom10.pdf")
