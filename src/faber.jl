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

function Faber(i, j)
    return Int(2^i + j)
end
function Faber_(n)
    return (Int(floor(log2(n))), Int(n - 2^(floor(log2(n)))))
end

const FaberArray = [Faber_(i) for i in 1:5000]
#You can precomile also Faber(i,j) as line above

Faber(n) = FaberArray[n]
