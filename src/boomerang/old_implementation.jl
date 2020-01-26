include("../ZZDiffusionBridge.jl")
using LinearAlgebra

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
    Q(α, ϕ, L, T)

generate matrix `Q` for the sin process X with diffusion length T.
Does not depend on ξ --> pre-compile
"""
function Q(α::Float64, ϕ::Vector{Fs}, L::Int64, T::Float64)
    dim = 2^(L+1) -1
    Q = zeros(dim,dim)
    for i in 1:dim
        for j in 1:i
            if i == j
                Q[i,i] = 0.5*ϕ[i].supr*ϕ[i].δ*(2*α^2 + α)
            else
                if ShareDomain(i,j)==1
                    ub = min(ϕ[i].ub, ϕ[j].ub)
                    lb = max(ϕ[i].lb, ϕ[j].lb)
                    S_ij = ub - lb
                    Q[i,j] = Q[j,i] =  0.5*S_ij*ϕ[i].supr*ϕ[j].supr*(2*α^2 + α)
                else
                    Q[i,j] = Q[j,i] = 0
                end
            end
        end
    end
    return Q
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

"""
    ∇U_0_numeric(α::Float64, L::Int64, T::Float64, u::Float64, v::Float64)

This has to be changed because I can compute the integral analytically.
"""
function ∇U_0_numeric(α::Float64, L::Int64, T::Float64, u::Float64, v::Float64)
    ∇U_0 = zeros(2^(L+1)-1)
    for n in 1:length(∇U_0)
        i, j = Faber(n)
        dt = 0.001
        dt2 = dt/2
        for t in 0:dt:T-dt
            ∇U_0[n] += 0.5*dt*Λ(t + dt2, i, j, T)*(α^2*sin(2*(u + (v-u)*(t + dt2)/T))
                        - α*sin((u + (v-u)*(t + dt2)/T)))
        end
    end
    return ∇U_0
end


#Poisson rates
"""
    wait_linear(a,b,u)
obtaining waiting time for Inhomogeneous Poisson Process
with rate of the form λ(t) = (a + b*t)^+, `a`,`b` ∈ R, `u` random variable
### CHANGE NAME
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
function event_λref(c::Float64)
    return -log(rand())/c
end

"""
    event_λbar(a::Float64, b::Float64)
Draw waiting time from rate λ(t) = max(0, a + b t)
"""
function event_λbar(a::Float64, b::Float64)
    wait_linear(a, b, rand())
end


"""
    ∇U_tilde(ξ::Vector{Float64}, ϕ::Vector{Fs}, α::Float64, L::Int64, T::Float64, u::Float64, v::Float64)

accept reject time drwan from upper bound λbar relative to the coefficient `n`
of model `SinSDE` starting at `u` and ending at `v`
"""
function ∇U_tilde(ξ::Vector{Float64}, ϕ::Vector{Fs}, α::Float64, L::Int64, T::Float64, u::Float64, v::Float64)
    ∇U_tilde = zeros(length(ϕ))
    for n in 1:length(ϕ)
        #δ = ϕ[n].range*(X.α*X.α + X.α)*0.5 Not used
        t = MCintegration(ϕ[n])
        XX = fs_expansion(t, ξ, ϕ, u, v, L, T)
        ϕ_t = Λ(t, ϕ[n].i, ϕ[n].j, T)
        ∇U_tilde[n] = 0.5*ϕ[n].range*ϕ_t*(α*α*sin(2.0*(XX)) - α*sin(XX))
    end
    return ∇U_tilde
end

function acc_rej(∇U_tilde::Vector{Float64}, θ::Vector{Float64})
    max(0, ∇U_tilde*θ)/(a + b*τ)
end



"""
    R(x, v , ∇U_tilde)

contour reflections
"""
function R(ξ::Vector{Float64}, θ::Vector{Float64} , ∇U_tilde::Vector{Float64})
    return θ - 2*dot(∇U_tilde, θ)/dot(∇U_tilde, ∇U_tilde)*∇U_tilde
end


function boomerang_sampler_sin(T::Float64, L::Int64, u::Float64, v::Float64, clock::Float64)
    clearconsole()
    ξ = zeros(2^(L+1)-1)
    θ = randn(2^(L+1)-1)
    t = 0.0
    α = 0.5
    c = 0.001
    ϕ = generate(L, T)
    norm_∇U_0 = norm(∇U_0_numeric(α, L, T, u, v)) ## maybe do it numerically now in a fine grid, later do it analytically.
    norm_Q =  opnorm(Q(α, ϕ, L, T))
    #Q_try = Q(α, ϕ, L, T) #DEBUG
    #if Q_try != Q_try'
    #    error("Q not symmetric")
    #end  # does not depend on x, v so precompile till here
    ∇Ubar = ∇U_bar(α, ϕ, L, T)# does not depend on x, v so precompile
    a = max(dot(abs.(θ), ∇Ubar),0)
    b = norm_∇U_0*sqrt(dot(θ, θ) + dot(ξ, ξ)) + norm_Q*(dot(θ, θ) + dot(ξ, ξ))
    τ0 = event_λbar(a, b)
    #println("time: ", τ0) #DEBUG
    #println("a + b*τ0: ", a + b*τ0) #DEBUG
    τ_ref = event_λref(c)
    Ξ = [Skeleton2(copy(ξ), copy(θ), t)]
    while t < clock
        if τ_ref < τ0
            #println("STEP: refreshment")
            ξ, θ =  boomerang_traj(ξ, θ, τ_ref)
            t += τ_ref
            θ = randn(length(θ))
            push!(Ξ, (Skeleton2(copy(ξ), copy(θ), t)))
            a = max(dot(abs.(θ), ∇Ubar),0)
            b = norm_∇U_0*sqrt(dot(θ, θ) + dot(ξ, ξ)) + norm_Q*(dot(θ, θ) + dot(ξ, ξ))
            τ0 = event_λbar(a, b)
            τ_ref = event_λref(c)
        else
            ξ, θ =  boomerang_traj(ξ, θ, τ0)
            t +=  τ0
            τ_ref -= τ0
            ∇Utilde = ∇U_tilde(ξ, ϕ, α, L, T, u, v)
            acc_ratio = max(dot(∇Utilde, θ), 0)/(a + b*τ0)
            if   !(0 <= acc_ratio < 1)
                #println("a + b*τ0: ", a + b*τ0)
                #println("lambda tilde ",  max(dot(∇Utilde, θ), 0))
                #println("acc_ratio: ", acc_ratio)
                error("invalid acc/rej ratio")
            end
            if acc_ratio > rand()
                #println("STEP: ok acc_ratio, accept event time")
                θ = R(ξ, θ, ∇Utilde)
                push!(Ξ, (Skeleton2(copy(ξ), copy(θ), t)))
                a = max(dot(abs.(θ), ∇Ubar),0)
                b = norm_∇U_0*sqrt(dot(θ, θ) + dot(ξ, ξ)) + norm_Q*(dot(θ, θ) + dot(ξ, ξ))
                τ0 = event_λbar(a, b)
            else
                #println("STEP: ok acc_ratio, reject event time")
                a = max(dot(abs.(θ), ∇Ubar),0)
                b = norm_∇U_0*sqrt(dot(θ, θ) + dot(ξ, ξ)) + norm_Q*(dot(θ, θ) + dot(ξ, ξ))
                τ0 = event_λbar(a, b)
            end
        end
    end
    return Ξ
end

T = 50.0
u = Float64(-π)
v =  Float64(3π)
L = 6
clock = 10000.0
Ξ = boomerang_sampler_sin(T, L, u, v, clock)

error("STOP HERE")

save("ouput.jld", "output", Ξ)


Ξ = load("ouput.jld")
Ξ = Ξ["output"]

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
    p = plot(dt, P, color = :darkrainbow, line_z = (1:N)', linewidth=0.001, alpha = 0.5, leg = false, colorbar = true)
    hline!(p, [n*π for n in -1:2:3], color = :blue)
    display(p)
    return p
end


using Plots
b = 20:10.0:9990
length(b)
plot_boomerang(Ξ, b, L, T, u, v)


savefig("sin_boomerang3.png")
