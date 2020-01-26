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

function plotmixing(y::Array{Skeleton,1}, b, T::Float64, L::Int64, u::Float64, v::Float64, hist, trasform = x -> x)
    # if hist
    #     ind = fill(0,L)
    #     [ind[i] =Faber(i,i) for i in 1:L]
    #     #dt = range(0.0, stop=T, length=2<<(L) + 1)
    #     P = []
    #     hist = zeros(length(b), length(ind))
    #     for i in b
    #         ξ_interp = FindCoordinates(y, i).ξ
    #         hist[j,:] = ξ_interp[ind]
    #         j = j + 1
    #     end
    #     return hist
    # else
    p = Plots.plot()
    P = []
    dt = range(0.0, stop=T, length=2<<(L) + 1)
    N = length(b)
    j = 1
    for i in b
        ξ_interp = FindCoordinates(y, i).ξ
        dx = fs_expansion(ξ_interp, u, v, L, T)
         push!(P, trasform.(dx))
    end
    Plots.plot!(p, dt, P, color = :darkrainbow, line_z = (1:N)' , linewidth=0.1, alpha = 0.2, leg = false, colorbar = true)
    return p
end



# function plotmixing(y::Array{Skeleton,1}, b, T::Float64, L::Int64, u::Float64, v::Float64, hist, trasform = x -> x)
#     if hist
#         ind = fill(0,L)
#         [ind[i] =Faber(i,i) for i in 1:L]
#         #dt = range(0.0, stop=T, length=2<<(L) + 1)
#         P = []
#         hist = zeros(length(b), length(ind))
#         for i in b
#             ξ_interp = FindCoordinates(y, i).ξ
#             hist[j,:] = ξ_interp[ind]
#             j = j + 1
#         end
#         return hist
#     else
#         p = Plots.plot()
#         dt = range(0.0, stop=T, length=2<<(L) + 1)
#         j = 1
#         for i in b
#             ξ_interp = FindCoordinates(y, i).ξ
#             dx = fs_expansion(ξ_interp, u, v, L, T)
#              push!(P, trasform.(dx))
#         end
#         Plots.plot(dt, P,  line_z = j, linewidth=0.3, alpha = 0.2, leg = false, colorbar = true)
#         return p
#     end
# end
