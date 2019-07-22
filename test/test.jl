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





function Λ(x, l, k, T)
   x = mod(x, T)   #ok
   2^(-l/2)*Λ((x)*1<<(l) -T*(k),T)
end

using Plots

T=1
dt= 0:0.0001:T
plot(dt, [Λ(t, 2, 0, T) for t in dt])





function num_int(i1,j1,i2,j2, T, dt)
    t = 0:dt:T
    r = 0
    for i in t[2:end]
        r += Λ(i + dt/2, i1, j1, T)*Λ(i + dt/2, i2, j2, T)*dt
    end
    return r
end


function num_int(i1,j1, T, dt)
    t = 0:dt:T
    r = 0
    for i in t[2:end]
        r += Λ(i + dt/2, i1, j1, T)*dt
    end
    return r
end



function num_int_bound1(n, T, dt)
    i1,j1 = Faber(n)
    t = 0:dt:T
    r = 0.0
    for i in t[2:end]
        r += Λ(i + dt/2, i1, j1, T)*(i + dt/2)/T*dt
    end
    return r
end


function num_int_bound2(n, T, dt)
    i1,j1 = Faber(n)
    t = 0:dt:T
    r = 0.0
    for i in t[2:end]
        r += Λ(i + dt/2, i1, j1, T)*(1 -(i + dt/2)/T)*dt
    end
    return r
end

num_int(0, 0, 2, 2 , 100, 0.00001)
num_int(0, 0 , 10, 0.00001)


generate_matrix(2, 10)
generate_vector(2, 10)


Λ(2.5, 1, 0, 1)
sqrt(10)*2^-(1/2 + 1)


bound2(3, 100)
num_int_bound2(3, 100.0, 0.001)
