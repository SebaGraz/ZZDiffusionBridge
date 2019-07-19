using Test

function Λ1(x, T)
    return (T-x)/T
end

#add basis for the final value
function Λ2(x,T)
    return x/T
end

function Λ(x, T)
    return max(0.0, sqrt(T)/2 - abs(x/sqrt(T) - sqrt(T)/2))
end
#T = 3
#plot([Λ(x, T) for x=0:0.001:T])

function Λ(x, l, k, T)
    x = mod(x, T)   #ok
    2^(-l/2)*Λ((x)*(2)^(l) - k*T, T)
end

function c(k::Fs, ξ, t, u, v, T)
    c = 0.0
    for n in k.nhb
        i,j = Faber(n)
        c += Λ(t, i, j, T)*ξ[n]
    end
    return c + Λ1(t, T)*u + Λ2(t, T)*v
end


L = 1
T = 10.0
ϕ = Fs(1, L, T)
ξ = fill(1.0, 2<<L - 1)
u = -1.0
v = -1.0




S = System(L, T, ξ)


fs_expansion(S, t, u, v)
c(ϕ, ξ, t, u, v, T)
for k in S.ϕ
    @test MCintegration(k) < k.ub
    @test MCintegration(k) > k.lb
    @test k.lb + range = k.ub
end

fs_expansion()
fs_expansion(ϕ, ξ, u, v, L, T)


for i
function MCintegration(k::Fs)
end









using Test
for i in 1:2^20
    @test i >= 1 << binlog(i)
    @test i < 2 << binlog(i)
    @test binlog(i) == Faber(i)[1]
end
