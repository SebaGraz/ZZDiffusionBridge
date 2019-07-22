"""
    OUSDE <: AbstractModel

dX_t = ν(μ - X_t)dt + dB_t
ν := intensity
μ := mean reversion
"""
struct OUSDE <: AbstractModel
    μ::Float64
    ν::Float64
end

dependence_strucute(::OUSDE) = PartialIndependence()
sampling_scheme(::OUSDE) = Regular()



#Poisson rates
"""
obtaining waiting time for Inhomogeneous Poisson Process
with rate of the form λ(t) = (a + b*t)^+, `a`,`b` ∈ R, `u` random variable
"""
function wait_gengaus(a,b,u)
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
    λbar(n, S::System, X::OUSDE , u, v)

Poisson time (real bound) for the coefficient `n`
of model `OUSDE` starting at `u` and ending at `v`

λbar = λ No upper bound and no accept reject step
"""
function λbar(n, S::System, X::OUSDE , u::Float64, v::Float64, t::Float64)
    a = S.θ[n]*(X.ν*X.ν*(dot(S.M[n,:], S.ξ) + S.bound1[n]*v + S.bound2[n]*u - X.μ*S.V[n]) + S.ξ[n])
    b = S.θ[n]*((dot(S.M[:,n], S.θ))*X.ν*X.ν + S.θ[n])
    return wait_gengaus(a, b, rand())
end

####

function runall(SHORT = false)
    T = 5.0
    clock = 100.0
    L = 5
    ν = 1.0
    μ = -5.0
    u = -1.0
    v = 1.0
    X = OUSDE(μ, ν)
    XX = zz_sampler(X, T, L, u, v, clock)
    if SHORT == false
        burning = 10.0    #burning
        f = clock - 1.0; n = 50
        db = (f-burning)/n
        b =  burning:db:f
        p = plotmixing(XX, b, T, L, u, v)
        display(p)
        #plot the mean of the process
    return XX
end

runall()
