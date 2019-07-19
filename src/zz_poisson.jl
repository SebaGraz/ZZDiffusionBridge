"""
    MCintegration(k::Fs)

take a random point inside the support of the basis function k
"""
function MCintegration(k::Fs)
    rand()*k.range + k.lb
end

"""
    Sol2E(a::Float64, b::Float64, c::Float64)

Solve a quadratic equation.
TODO add errors when the solution is not real.
"""
#solve quadratic equation ax^2 + bx + c = 0
function Sol2E(a::Float64, b::Float64, c::Float64)
    return (-b + sqrt(b^2 - 4a*c))/2a
end

######################################
##############  Sin sde     ##########
######################################
"""
        λbar(n, S::System, X::SinSDE, u, v)

Poisson time (upper bound) for the coefficient `n`
of model `SinSDE` starting at `u` and ending at `v`

invert function Λ(t) = at + int_0^t max(0, b+cs) ds + ln(ran)
with a = |θ|*δ,   b = ξ*θ, c = θ^2
"""
function λbar(n, S::System, X::SinSDE, u, v, t::Float64)
    ran = rand()
    δ = S.ϕ[n].δ*(X.α*X.α + X.α)*0.5   #always the same, we could save the value
    b = S.ξ[n]*S.θ[n]
    if b>0
        return Sol2E(S.θ[n]*S.θ[n]*0.5, b + δ*abs(S.θ[n]), log(ran))
    elseif S.ξ[n]*sign(S.θ[n])*δ <= log(ran) #Case 2: 0 < t < -b/c
        return -log(ran)/(abs(S.θ[n])*δ)
    else    #Case 2: 0 < -b/c < t
        return Sol2E(S.θ[n]*S.θ[n]*0.5, b+δ*abs(S.θ[n]), log(ran)+ S.ξ[n]*S.ξ[n]*0.5)
    end
end

"""
    λratio(n::Int64, S::System, X::SinSDE, u::Float64, v::Float64)

accept reject time drwan from upper bound λbar relative to the coefficient `n`
of model `SinSDE` starting at `u` and ending at `v`
"""
function λratio(n::Int64, S::System, X::SinSDE, u::Float64, v::Float64, t::Float64)
    δ = S.ϕ[n].δ*(X.α*X.α + X.α)*0.5 #always the same, we could save the value
    t = MCintegration(S.ϕ[n])
    XX = fs_expansion(S, t, u, v)
    λ = Λ(t, S.ϕ[n].i , S.ϕ[n].j, S.T)
    return max(0, S.θ[n]*(0.5*S.ϕ[n].range*λ*(X.α*X.α*sin(2.0*(XX)) - X.α*sin(XX)) + S.ξ[n]))/(abs(S.θ[n])*δ + max(0, S.θ[n]*S.ξ[n]))
end

###################################
####    LOGISTIC GROWTH SDE  ######
###################################
"""
    waitgaus(ξ, θ, u)

Poisson rate  λ(s) = max(0, θ(ξ + θs))
invert function Λ(t) = int_0^t λ(s) ds + ln(u)
"""
function waitgaus( ξ::Float64, θ::Float64, u::Float64)
    b = ξ*θ
    if b>0
        return Sol2E(0.5, b, log(u))
    else
        return Sol2E(0.5, b, log(u) + ξ*ξ*0.5)
    end
end
"""
    waitexp(a, b, u)

invert function  Λ(t) = a/b (e^(bt) - 1) + ln(u)
"""
function waitexp(a, b, u)
    if (log(u)*b/a > 1)
        return error("exponential gave infinity")
    else
        return log(1-log(u)b/a)/b
    end
end
"""
    λbar(n, S::System, X::LogGrowthSDE , u, v)

Poisson time (upper bound) for the coefficient `n`
of model `LogGrowthSDE` starting at `u` and ending at `v`
"""
function λbar(n, S::System, X::LogGrowthSDE , u::Float64, v::Float64, t::Float64)
    w1 = waitgaus(S.ξ[n], S.θ[n], rand())
    S.b1[n] = minimum(fs_expansion(S.ϕ[n], S.ξ, u, v, S.L, S.T))
    S.b2[n] = minimum(fs_expansion(S.ϕ[n], S.θ, u, v, S.L, S.T))
    S.tt[n] = t
    if S.θ[n] > 0
        w2 = waitexp(0.5*S.ϕ[n].δ*X.a1*exp(-X.β*S.b1[n]), -X.β*S.b2[n], rand())
        return min(w1,w2)
    else
        w2 = waitexp(0.5*S.ϕ[n].δ*X.a2*exp(-2X.β*S.b1[n]), -2X.β*S.b2[n], rand())
        return min(w1, w2)
    end
end


"""
    λratio(n::Int64, S::System, X::LogGrowthSDE, u::Float64, v::Float64)

accept reject time drwan from upper bound λbar relative to the coefficient `n`
of model `LogGrowthSDE` starting at `u` and ending at `v`
"""
function λratio(n::Int64, S::System, X::LogGrowthSDE, u::Float64, v::Float64, t::Float64)
    time = t - S.tt[n]
    ω = MCintegration(S.ϕ[n])
    XX = fs_expansion(S, ω, u, v)
    ϕ = Λ(ω, S.ϕ[n].i, S.ϕ[n].j, S.T)   #need to be changed
    num = max(0, S.θ[n]*(0.5*S.ϕ[n].range*ϕ*(X.a1*exp(-X.β*XX) - X.a2*exp(-2X.β*XX)) + S.ξ[n]))
    den = max(0, S.θ[n]*S.ξ[n]) + max(0, 0.5*S.θ[n]*S.ϕ[n].δ*X.a1*exp(-X.β*S.b1[n])*exp(-X.β*S.b2[n]*time)) + max(0, -0.5*S.θ[n]*S.ϕ[n].δ*X.a2*exp(-2X.β*S.b1[n])*exp(-2X.β*S.b2[n]*time))
    return num/den
end





########################
####    OU Model    ####
########################

"""
    λbar(n, S::System, X::OUSDE , u, v)

Poisson time (real bound) for the coefficient `n`
of model `OUSDE` starting at `u` and ending at `v`

λbar = λ No upper bound and no accept reject step
"""
function λbar(n, S::System, X::OUSDE , u, v)

end
