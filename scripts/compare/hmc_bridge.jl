##### Packaged loglikelihood
using LinearAlgebra
#using Turing
using Random
using AdvancedHMC
using ForwardDiff
include("../../src/faberschauder.jl")
include("../../src/asvar.jl")
include("../../src/pdmp.jl")
# up to constant of proportionality
# π(x) = exp(-ϕ)*C
function ϕ(ξ, L, T; p = T/(1<<L)) # formula (17)
        r = 0.0
        s = p*0.5
        for i in 1:(1<<L)
            x = dotψ(ξ, s, L,  T)
            r += 0.5*p*(b(x)^2 + b′(x))
            s += p
        end
        @assert  T+p*0.5 - 0.0001 < s < T+p*0.5 + 0.0001
    return r + 0.5*norm(ξ)^2
end

function ∇ϕ(ξ, i, L, T; p = T/(1<<L)) # formula (17)
    if i == (2 << L) + 1    # final point
        error("fixed final point")
    elseif i == 1   # initial point
        error("fixed initial point")
    else
        l = lvl(i, L)
        k = (i - 1) ÷ (2 << l)
        # println("index ($l, $k) out of $L level")
        δ = T/(1 << (L-l)) # T/(2^(L-l))
        r = 0.0
        s = δ*k + p*0.5
        # println("T = $T")
        for i in 1:l+1
            # println("integrating at the point $s")
            x = dotψ(ξ, s, L,  T)
            r += 0.5*p*Λ(s, L-l, T)*(2b(x)*b′(x) + b″(x))
            s += p
        end
        @assert δ*(k + 1) + p*0.5 - 0.001 < s < δ*(k + 1) + p*0.5 + 0.001
        return r + ξ[i]
    end
end


# Drift
b(x) = α * sin(x)
# First derivative
b′(x) = α * cos(x)
# Second derivative
b″(x) = -α * sin(x)

# Drift
#α = 1.5
L = 6
T = 50.0
n = (2 << L) + 1
u, v = 0.0, 0.0  # initial and fianl point
D = n - 2

# loglikelihood
function logdensity_f(q)
    ξ = [u/sqrt(T); q;  v/sqrt(T)]
    return -ϕ(ξ, L, T)
end
# gradient of loglikelihood
function grad_f(q)
    ξ = [u / sqrt(T); q; v / sqrt(T)]
    return (-ϕ(ξ, L, T),[-∇ϕ(ξ, i, L, T) for i in 2:(2 << L)])
end

# reshape output from vector of vector to matrix
function vecvec_to_matrix(vecvec)
    dim1 = length(vecvec)
    dim2 = length(vecvec[1])
    my_array = zeros(Float64, dim2, dim1)
    for i in 1:dim1
        for j in 1:dim2
            my_array[j,i] = vecvec[i][j]
        end
    end
    return my_array
end



function run_hmc(df)
    q0 = randn(D)
    n_samples, n_adapts, target = 3_000, 2_000, 0.8
    metric = DiagEuclideanMetric(D)
    # h = Hamiltonian(metric, logdensity_f, grad_f)
    h = Hamiltonian(metric, logdensity_f, ForwardDiff)
    eps_init = find_good_stepsize(h, q0)
    int = Leapfrog(eps_init)
    traj = AdvancedHMC.NUTS{MultinomialTS,GeneralisedNoUTurn}(int)
    adaptor = AdvancedHMC.StanHMCAdaptor(
        n_adapts, MassMatrixAdaptor(metric), NesterovDualAveraging(target, eps_init)
        )
    #u, v = 0.0/sqrt(T), 0.0/ sqrt(T)  # initial and fianl point
    hmc_time = @elapsed ((samples, stats) = sample(h, traj, q0, n_samples, adaptor, n_adapts))
    burn1 = 100
     # xT_2 = sqrt(T)/2*(samples[Int((length(q0)+1)/2), burn1+1:niter])
    xT_2 = sqrt(T)/2 .*(getindex.(samples, Int((length(q0)+1)/2)))
    nn = length(xT_2)
    ave = cumsum(xT_2)./(Float64.(1:nn))
    for batches in [50]
        #reshape a vector of vectors in a matrix
        ess = ESS(vecvec_to_matrix(samples)[:,burn1:end], n_batches = batches)
        push!(df, Dict(:sampler => "HMC", :alpha => α, :T => n_samples, :nbatches => batches,
            :stat => "ess_mean", :y => sum(ess[2:end-1])/(length(ess)-2), :runtime =>  hmc_time))
        push!(df, Dict(:sampler => "HMC", :alpha => α, :T => n_samples, :nbatches => batches,
            :stat => "ess_median", :y => median(ess[2:end-1]), :runtime =>  hmc_time))
        push!(df, Dict(:sampler => "HMC", :alpha => α, :T => n_samples, :nbatches => batches,
            :stat => "ess_min", :y => minimum(ess[2:end-1]), :runtime => hmc_time))
        push!(df, Dict(:sampler => "HMC", :alpha => α, :T => n_samples, :nbatches => batches,
            :stat => "ess_x_T2", :y => ess[Int((length(q0)+1)/2)], :runtime => hmc_time))
    end
    return df, ave
end


using DataFrames, Random
function data_collection_hmc(df)
    Random.seed!(0)
    for α′ in [0.1, 0.3, 0.5]
        global α = α′
        println("iteration ")
        df, x = run_hmc(df)
    end
    return df
end

df = DataFrame(sampler = String[], alpha = Float64[], T = Float64[],
    nbatches = Int64[], stat = String[], y = Float64[], runtime = Float64[])
df = data_collection_hmc(df)
using CSV, JLD2
CSV.write("./scripts/compare/benchamrk_hmc.csv", df)
