using Plots
using LaTeXStrings

function interpolate(x1, y1, x2, y2, x)
    return ( y1 .+ (y2 .- y1)/(x2 .- x1).*(x .- x1))
end


function FindCoordinates(Ξ::Array{Skeleton,1}, t)
    tt = [Ξ[i].t for i in 1:length(Ξ)]
    i = searchsortedfirst(tt, t)
    #i = findfirst(x -> x>t ,tt)
    return(Skeleton(interpolate(Ξ[i-1].t, Ξ[i-1].ξ, Ξ[i].t, Ξ[i].ξ, t), t))
end



#x1 = Skeleton(ones(2<<1- 1), 0.0)
#x2 = Skeleton(fill(0, 2<<1 - 1), 1.0)
#FindCoordinates([x1,x2], 0.5).ξ

function plotmixing(y::Array{Skeleton,1}, b, T::Float64, L::Int64, u::Float64, v::Float64, trasform = x -> x)
    p = Plots.plot(leg = false, colorbar = true)
    dt = range(0.0, stop=T, length=2<<(L) + 1)
    j = 1
    for i in b
        dx = fs_expansion(FindCoordinates(y, i).ξ, u, v, L, T)
        Plots.plot!(p, dt, trasform.(dx),  line_z = j, linewidth=0.3, alpha = 0.2)
        j = j + 1
    end
    return p
end
