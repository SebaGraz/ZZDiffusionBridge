
using Plots
get_info(x) = x, norm(x)^2/2




function mala_sampler(niter::Int64, ξ_init =fill(0.0, 2))
        ξ₀ = ξ_init
        ξtrace=Array{Float64}(undef, length(ξ_init), niter)
        ξtrace[:,1] = ξ₀
        ∇U₀, U₀ = get_info(ξ₀)
        τ = 1.0
        count = 0
        for i=2:niter
                ξ₁  = ξ₀ .- τ*∇U₀ .+ sqrt(2*τ)*randn(length(ξ₀))
                ∇U₁, U₁ = get_info(ξ₁)
                acc_rej =  U₀ - U₁ + (norm(ξ₁ - ξ₀ + τ*∇U₀)^2 - norm(ξ₀ - ξ₁ + τ*∇U₁)^2)/(4τ)
                if acc_rej > log(rand())
                        ∇U₀, U₀, ξ₀ = ∇U₁, U₁, ξ₁
                        count += 1
                end
                # daptive step here
                ξtrace[:,i] = ξ₀
        end
        println("print tot ar : $(count/niter)")
        ξtrace
end


y = mala_sampler(50000)

scatter(y[1,:],y[2,:], aspect = :equal, alpha = 0.01)
plot!([(2*sin(x), 2*cos(x)) for x in -π:0.001:pi], color = :red, label = "", linewidth = 2, )
using Statistics
mean(y[1,:])


50/(2<<6)
