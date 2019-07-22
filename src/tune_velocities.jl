using BenchmarkTools
include("ZZDiffusionBridge.jl")

#comment out runall() in example_sin.jl
include("./examples/example_sin.jl")

# see notes section: tune velocities
function tune_velocity( L::Int64, β = 0.0)
    θ = fill(0.0, 2<<L - 1)
    k = 1
    for i in 0:L
        for j in 0:2^i-1
            θ[k] = 2.0^(-i*β)
            k += 1
        end
    end
    return(θ)
end


#just little differences from the original zz_sampler (src/zz_sampler).
function zz_sampler(S::System, X::AbstractModel, T::Float64, L::Int64, u::Float64, v::Float64, clock::Float64, τ0::Float64 , n0::Int64)
    t = 0.0
    Ξ = Skeleton[]
    while t < clock
        τ0 , n0 = findmin(S.τ)
        S.ξ .+=  S.θ*τ0
        t += τ0
        if acc_rej(n0, S, X, sampling_scheme(X), u, v, t)
            S.θ[n0] = -S.θ[n0]
            push!(Ξ, (Skeleton(copy(S.ξ), t)))
            update_events!(n0, τ0, S, X, sampling_scheme(X), dependence_strucute(X), true, u, v, t)
        else
            update_events!(n0, τ0, S, X, sampling_scheme(X), dependence_strucute(X), false, u, v, t)
        end
    end
    Ξ, τ0, n0
end

function piec_lin_integration(y, L)
    N = 2<<L
    W = 1/N
    somma = 0.5*(y[1] + y[end])*W + sum(y[2:(end-1)])*W
end


function compute_average(X::Array{Skeleton,1}, sampling_frequncy::Float64, u, v, L)
    T = X[end].t - 1.0 #0.1 prevents error in the interpolation
    dt = 1.0:sampling_frequncy:T
    M = 0
    count = 0
    for i in dt
        ξ = FindCoordinates(X, i).ξ
        ξfe = fs_expansion(ξ, u, v, L, T)
        path = piec_lin_integration(ξfe, L)
        M = M + path
        count += 1
    end
    return M/count
end

"""


"""
function vel_experiment(β, L, clock)
    α = 0.0
    T = 200.0
    u = -3π
    v = 3π
    τ0 = 0.0
    n0 = 0
    ξ = fill(0.0, 2<<L - 1)
    θ = tune_velocity(L, β)
    S = System(L, T, ξ, θ)
    X = SinSDE(α)
    Y = []
    TIME = []
    tot_time = 0
    for i in 1:20
        (Ξ, τ0, n0), time = @timed  zz_sampler(S, X, T, L, u, v, clock, τ0, n0)
        #S::System, X::AbstractModel, T::Float64, L::Int64, u, v, clock, τ0 , n0)
        tot_time += time
        push!(Y, copy(Ξ) )
        push!(TIME, tot_time )
    end
    return Y, TIME
end

function piec_lin_integration(y, L)
    N = 2<<L
    W = 1/N
    somma = 0.5*(y[1] + y[end])*W + sum(y[2:(end-1)])*W
end




function runall_velocity()
    Y = fill(0.0,(6,20))
    Time = fill(0.0,(6,20))
    row = 1
    L = 7
    p = plot()
    for β in 0.0:0.1:0.5
        X, time = vel_experiment(β, L, 500.0)
        ave = cumsum(compute_average.(X, 1.0, -3π,  3π, L))
        Y[row,:] =  [ave[i]/i for i in 1:20]
        Time[row, :] = time
        plot!(p, Time[row,:], abs.(Y[row,:]), label =("beta =  $β"))
        row += 1
    end
    display(p)
    return Y, Time
end



X, T = runall_velocity()
