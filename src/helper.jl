#∑ ϕ_i ξ_i for i sharing the same support of the basis k
function c(k::GFSbase, ξ, t, u, v)
    c = 0.0
    for n in k.nhb
        i,j = Faber(n)
        c += Λ(t, i, j, k.T)*ξ[n]
    end
    return c + Λ1(t, k.T)*u + Λ2(t, k.T)*v
end

#sample randomly from the support of the basis
function MCintegration(k::GFSbase, n = 1)
    rand(n).*k.range .+ k.lb
end

#solve quadratic equation ax^2 + bx + c = 0
function Sol2E(a::Float64, b::Float64, c::Float64)
    return (-b + sqrt(b^2 - 4a*c))/2a
end

#invert equation ???? Λ(t) = δ*θ*t + int_0^t |ξ + θ*s| ds + ln(u)
function TimeAbs(k::GFSbase, ξ, θ, α, σ, u)
    δ = k.δ*α*(α/σ + 1)*0.5
    if ξ*θ >= 0
        θ = abs(θ)
        ξ = abs(ξ)
        return Sol2E(θ^2/2, θ*(δ + ξ), log(u))
    else
        ξ = - abs(ξ)
        θ = abs(θ)
        #if  -δ*ξ + 0.5*ξ^2 <= -log(u)
        if  -δ*ξ + 0.5*ξ^2 <= -log(u)
            return Sol2E(θ^2/2, θ*(δ + ξ), ξ^2 + log(u))
        else
            return Sol2E(-θ^2/2, θ*(δ - ξ), log(u))
        end
    end
end

function λratio(k::GFSbase, ξ, θ, τ, α, σ, u, v, N) #sta attento a non mpdificare il vettore
    ϵ = 0
    t = MCintegration(k, N)
    w = 1/N
    for i in 1:N
        cc = c(k, ξ, t[i], u, v)
        ϕ = Λ(t[i], k.i ,k.j, k.T)
        ϵ += w*ϕ*(α/σ*sin(2.0*(cc + ξ[k.n]*ϕ )) - sin(cc + ξ[k.n]*ϕ ))
    end
    return max(0, θ[k.n]*(0.5*k.range*α/σ*ϵ + ξ[k.n]))/((α/σ + 1)*0.5*k.δ*α*abs(θ[k.n]) + abs(ξ[k.n]*θ[k.n]))
end
