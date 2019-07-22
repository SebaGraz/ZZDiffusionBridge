include("faber.jl")
include("types.jl")


"""
    acc_rej(n::Int64, S::System, X::AbstractModel, f1::Regular, u, v)

For Regular sampling_scheme: always accept
"""
function acc_rej(n::Int64, S::System, X::AbstractModel, f1::Regular, u::Float64, v::Float64, t::Float64)
    return true
end
"""
    acc_rej(n::Int64, S::System, X::AbstractModel, f1::SubSampling, u, v)

Subsampling scheme: accept the event time as coming from the real Poisson rate.
"""
function acc_rej(n::Int64, S::System, X::AbstractModel, ::SubSampling, u::Float64, v::Float64, t::Float64)
    if λratio(n, S, X, u, v, t)> rand()
        return true
    else
        return false
    end
end

"""
    update_events!(n, S::System, X::AbstractModel, ::SubSampling, ::FullIndependence, acc)

renovate time when the nth coefficient rings (i.e. τ_n = 0). In the full independece case
we do not need to simulate any other time, but just rescale them.
"""
function update_events!(n::Int64, τ0::Float64, S::System, X::AbstractModel, ::SubSampling, ::FullIndependence, acc::Bool, u::Float64, v::Float64, t::Float64)
    S.τ[n] = λbar(n, S, X , u, v, t)
    for i in 1:(n-1)
        S.τ[i] -= τ0
    end
    for i in  (n + 1): length(S.ϕ)
        S.τ[i] -=  τ0
    end
end


"""
    update_events!(n, S::System, X::AbstractModel, f1::Regular, f2::PartialIndependence, acc)

renovate time when the `n`_th coefficient rings (i.e. τ_n = 0). In the partial independece case
we need to simulate all the time relative to the functions sharing the same support and rescale the others.
"""
function update_events!(n::Int64, τ0::Float64, S::System, X::AbstractModel, f1::Regular, f2::PartialIndependence, acc::Bool, u::Float64, v::Float64, t::Float64)
    for i in S.ϕ[n].nhb
        S.τ[i] = λbar(i, S, X , u, v, t)
    end
    for i in S.ϕ[n].notnhb
        S.τ[i] -= τ0
    end
end

"""
    update_events!(n, S::System, X::AbstractModel, f1::SubSampling, f2::PartialIndependence, acc)

renovate time when the `n`_th coefficient rings (i.e. τ_n = 0). In the partial independece case
we need to simulate all the time relative to the functions sharing the same support and rescale the others.
if `acc` is false, means that in the previous step we rejected the event and we did not change velocity.
This implies that we need to renovate only the `n`th time like in the Full independece case.
"""
function update_events!(n::Int64, τ0::Float64, S::System, X::AbstractModel, ::SubSampling, ::PartialIndependence, acc::Bool, u::Float64, v::Float64, t::Float64)
    if acc == true
        update_events!(n, τ0, S, X, Regular(), PartialIndependence(), acc, u, v, t)
    else
        update_events!(n, τ0, S, X, SubSampling(), FullIndependence(), acc, u, v, t)
    end
end



"""
    zigzagsampler(X::AbstractModel, T, L, u, v, TT)

run the ZigZag sampler for diffusion `X` conditioned to start at `u` and
end at `v` at time `T`. The infinite summation is truncated at level `L` and
the ZigZag process will run up to time `TT`. Optionally set velocities θ
"""
function zz_sampler(X::AbstractModel, T::Float64, L::Int64, u::Float64, v::Float64, clock::Float64; θ = fill(1.0, 2<<L - 1))
    t = 0.0
    S = System(L, T)
    τ0 = 0.0
    n0 = 0
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
    return Ξ
end
