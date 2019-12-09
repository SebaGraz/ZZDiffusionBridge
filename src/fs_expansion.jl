
"""
        fs_expansion(x, u, v, L, T, n = i -> 2^-(1 + i/2))
Global recursion: find the value of the process on the dyadic points given by the first L Levels,
`x` is a vector containing coordinates `u`,`v` initial and final point
`T` length of the basis
`n` is the normalization function (in our case known)
"""
function fs_expansion(x::Vector{Float64}, u::Float64, v::Float64, L::Int64, T::Float64, n = i -> 2^-(1 + i/2))
        y = fill(NaN, 2<<L + 1)
        y[1] = u
        y[end] = v
        k = 1
        c = sqrt(T)
        for i in 0:L
                for j in 0:(1<<i - 1)
                        l = 1<<(L-i) + 2<<(L-i)*(j) + 1 #new index
                        y[l] = x[k]*n(i)*c + (y[l + 1<<(L-i)] + y[l - 1<<(L-i)])/2
                        k += 1
                end
        end
        y
end




"""
        fs_expansion(ϕ::Fs, ξ, u, v, L::Int64, T::Float64, n = i -> 2^-(1 + i/2))
Local recursion, like Global but on the support of the Faber Schauder function `ϕ` up to level `L`
This function evaluates only the necessary functions for the evaluation of the points inside the support
of the function ϕ.
`x` is a vector containing coordinates `u`,`v` initial and final point
`T` length of the basis
`n` is the normalization function (in our case known)

"""
function fs_expansion(ϕ::Fs, ξ::Vector{Float64}, u::Float64, v::Float64, L::Int64, T::Float64, n = i -> 2^-(1 + i/2))
        y = fill(NaN, 2<<(L) + 1)   #number of points in which we require the evaluation of the function
        y[end] = v
        y[1] = u
        c = sqrt(T)
        for k in ϕ.nhb
                i,j = Faber(k)
                l = 2<<(L-i-1) + 2<<(L-i)*(j) + 1 #new index
                y[l] = ξ[k]*n(i)*c + (y[l + 2<<(L-i-1)]  + y[l - 2<<(L-i-1)])/2
                k += 1
        end
        y[(2<<(L-ϕ.i)*ϕ.j + 1):(2<<(L-ϕ.i)*(ϕ.j + 1) + 1)]
end




"""
        interpolate(x:: array, y::array, x0)

simple interpolation function given `x`, `y`: 2-element array, `x0` floating
"""
function interpolate(x::Vector{Float64}, y::Vector{Float64}, x0::Float64)
        if x0 == x[1]
                return y[1]
        elseif x0 == x[2]
                return y[2]
        else
                return ( y[2] - (y[2] - y[1])/(x[2] - x[1])*(x[2] - x0))
        end
end


"""
        fs_expansion(S::System, t::Float64, u, v, n = i -> 2^-(1 + i/2))

Find value of the process (Piecewise linear) for any t ∈ [0, T]
does not evaluate the  whole path, but just the points needed for t
"""
function fs_expansion(S::System, t::Float64, u::Float64, v::Float64, n = i -> 2^-(1 + i/2))
        dt = 0:S.T/(2<<S.L):S.T
        k = (searchsortedfirst(dt, t) - 1)
        j0 =  Int(ceil(k/2))-1
        n0 = Faber(S.L, j0)
        if k % 2 != 0
                return interpolate([dt[k], dt[k + 1]], fs_expansion(S.ϕ[n0], S.ξ, u, v, S.L, S.T)[1:2], t)
        else
                return interpolate([dt[k], dt[k + 1]], fs_expansion(S.ϕ[n0], S.ξ, u, v,  S.L, S.T)[2:3], t)
        end
end
