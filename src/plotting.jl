using Plots

function level(i)
    j = 0
    while (i != 0)
        i = i >> 1
        j = j + 1
    end
    return j-1
end

function interpolate(x1, y1, x2, y2, x)
    return ( y1 .+ (y2 .- y1)/(x2 .- x1).*(x .- x1))
end

function FindCoordinates(Ξ::Array{Skeleton,1}, t::Float64)
    tt = [Ξ[i].t for i in 1:length(Ξ)]
    i = searchsortedfirst(tt, t)
    #i = findfirst(x -> x>t ,tt)
    return(Skeleton(t, interpolate(Ξ[i-1].t, Ξ[i-1].ξ, Ξ[i].t, Ξ[i].ξ, t)))
end

#Plot Skeleton and choose the dimension
function Plots.plot(Ξ::Array{Skeleton,1},dims::Vector)
    x = [Ξ[i].ξ[dims[1]] for i in 1 : length(Ξ)]
    y = [Ξ[i].ξ[dims[2]] for i in 1 : length(Ξ)]
    plot(x,y)
end

#### Plots
function cies(t, ξ, L, T, u, v)
    x = 0.0
    i = 1
    for l in 0:L
        for k in 0:2^(l)-1
            x += ξ[i]*Λ(t, l, k, T)        #### Can be much more efficient, sum of function having support on t
            i += 1
        end
    end
    x + Λ1(t, T)*u + Λ2(t, T)*v
end


function plotmixing(y::Array{Skeleton,1}, b, T, L, u, v)
    p = plot()
    dt = range(0.0, stop=T, length=2^(L+1))
    for i in b
        dx = [cies(ti, FindCoordinates(y, i).ξ, L, T, u, v) for ti in dt]
        plot!(p, dt, dx,  color= RGB(i/b[end], 0.4, 1-i/b[end]), linewidth=0.3, alpha = 0.4)
    end
    plot!(p, leg=false)
    hline!(p, [n*π for n in -5:2:5])
    display(p)
end


function makegif(y::Array{Skeleton,1}, b, T, L, u, v)
    dt = range(0.0, stop=T, length=2^(L+1))
    k = 1
    for i in b
        dx = [cies(ti, FindCoordinates(y, i).ξ, L, T, u, v) for ti in dt]
        plot(dt, dx, linewidth=0.4, alpha = 0.6, leg = false)
        hline!([n*π for n in -5:2:5])
        savefig("images/$k")
        k += 1
    end
end
