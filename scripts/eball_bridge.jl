using Makie
using Random
#ϵ ball bridge in likely points with Euler discretization
include("./examples/sin.jl")
#drift(X::LogGrowthSDE, )
drift(X::SinSDE, x) = X.α*sin(x)

function ϵballBridge(S, ϵ, T,discretization, u, v, N)
    Int(T / discretization)
    Ξ = Vector{Float64}[]
    time_grid = 0:discretization:T
    x = fill(0.0, length(time_grid))
    x[1] = u
    while length(Ξ)<N
        for i in 1:length(time_grid)-1
            x[i+1] = x[i] +  drift(S, x[i])*discretization + sqrt(discretization)randn()
        end
        if x[end] < v + ϵ && x[end] > v - ϵ
            push!(Ξ, copy(x))
        end
    end
    return Ξ
end

Random.seed!(1)
ebridge = ϵballBridge(SinSDE(.7,4,5),0.05,50.0,0.0005,-π,π, 500)
#a = plot(0:0.0005:50, ebridge[1:30], leg=false)
#hline!(a, [n*π for n in -3:2:5])
a = Float64[]
push!(a, 1.0)
function scatter_joint(a, ebridge)
    x1 = Float64[]
    x2 = Float64[]
    for i in 1:length(ebridge)
        push!(x1, ebridge[i][Int((length(ebridge[i]) - 1)/4) + 1])
        push!(x2, ebridge[i][Int((length(ebridge[i]) - 1)*3/4) + 1])
    end
    Makie.scatter!(a, x1, x2, color = (:darkred, 1.0), markersize = 0.25)
end
using Plots
y = runall(true)

l1 = Λ(50/4, 0, 0, 50)
l2 = Λ(50*3/4, 1, 1, 50)

y1 = [1/4*π + 3/4*(-π) + y[i].ξ[1]*l1 + y[i].ξ[2]*l2 for i in 1:length(y)]
y2 = [1/4*(-π) + 3/4*π + y[i].ξ[1]*l1 + y[i].ξ[3]*l2 for i in 1:length(y)]


using Makie
ed = lines(y1, y2, color = (:lightblue), linewidth=0.5)
lines!(ed, y1, y2, color = (:blue, 0.01), linewidth=1.0)
scatter_joint(ed, ebridge )
axis = ed[Axis]
axis[:names, :axisnames] = ("X_12.5", "X_37.5")
Makie.save("joint_distribution.png", ed)

Makie.scatter!
