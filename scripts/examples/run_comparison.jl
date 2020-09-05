
include("./sin_zz.jl")
time_zz, XX_zz = runall_zz()
include("./sin_bouncy.jl")
time_bouncy, XX_bouncy = runall_bouncy()
include("./sin_mala.jl")
time_mala, XX_mala = runall_mala()
include("./sin_mala_pCN.jl")
time_pcn_mala, XX_pcn_mala = runall_pcn_mala()
@show time_zz, time_bouncy, time_mala, time_pcn_mala


function trajectory(y, crd)
    burning = 1.0
    db = 0.1
    f_clock = y[end].t - 0.01
    CC = Vector{Float64}()
    b =  burning:db:f_clock
    j = searchsortedfirst(getfield.(y, :t), burning) - 1
    for ti in b
        while ti > y[j+1].t
            j = j+ 1
        end
        ξᵢ = interpolate(y[j].t, y[j].ξ, y[j+1].t, y[j+1].ξ, ti)[crd]
        push!(CC, ξᵢ)
    end
    CC
end


p1 = plot(plot(trajectory(XX_zz, 1)), plot(trajectory(XX_zz, 2)),
        plot(trajectory(XX_zz, 3)), layout = (3,1))
savefig("output/zz.pdf")

p2 = plot(plot(trajectory(XX_bouncy, 1)), plot(trajectory(XX_bouncy, 2)),
        plot(trajectory(XX_bouncy, 3)), layout = (3,1))
savefig("output/bps.pdf")

p3 = plot(plot(XX_mala[1,:]), plot(XX_mala[2,:]), plot(XX_mala[3,:]),
        layout = (3,1))
savefig("output/mala.pdf")

p4 = plot(plot(XX_pcn_mala[1,:]), plot(XX_pcn_mala[2,:]),
        plot(XX_pcn_mala[3,:]), layout = (3,1))
savefig("output/mala_pcn.pdf")


p5 = plot(plot(trajectory(XX_zz, 1), label = "zz", ylim = (-4,4)), plot(trajectory(XX_bouncy, 1), ylim = (-4,4), label = "bps"),
            plot(XX_mala[1,:], ylim = (-4,4),label = "mala"), plot(XX_pcn_mala[1,:], ylim = (-4,4), label = "pcn_mala"), layout = (2,2))

savefig("output/comparison.pdf")
