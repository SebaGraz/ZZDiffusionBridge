# dependence structure for each stochastic differential equation
dependence_strucute(::LogGrowthSDE) = PartialIndependence()
dependence_strucute(::SinSDE) = FullIndependence()
dependence_strucute(::OUSDE) = PartialIndependence()


# sampling scheme for each stochastic differential equation
sampling_scheme(::LogGrowthSDE) = SubSampling()
sampling_scheme(::SinSDE) = SubSampling()
sampling_scheme(::OUSDE) = Regular()


function runall(sin = true)
    T = 200.0
    clock = 200.0
    L = 7
    if sin == true
        α = 0.7 #sin
        u = - 3π
        v = + 3π
        X = SinSDE(α)
    else
        K = 2000
        r = 0.1
        β = 0.1
        u = -log(50)/β     # end points in the lamperti transform
        v= -log(1000)/β
        X = LogGrowthSDE(r, K, β)   #end points in the lamperti tranform
    end
    XX = zigzagsampler(X::AbstractModel, T, L, u, v, clock)
    burning = 10.0    #burning
    f = clock - 1.0; n = 30
    db = (f-burning)/n
    b =  burning:db:f
    if sin == true
        p = plotmixing(XX, b, T, L, u, v)
        hline!(p, [n*π for n in -5:2:5])
    else
        p = plotmixing(XX, b, T, L, u, v, x -> exp(- x*β))
        hline!(p, [K])
    end
    display(p)
    return XX
end

x = runall(false)
error("STOP HERE")
