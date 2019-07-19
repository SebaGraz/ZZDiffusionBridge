"""
    Faber(i, j)

change of indexation: from double indexation i,j to single indexation n
"""
function Faber(i::Int64, j::Int64)
    return Int(2<<(i-1) + j)
end

"""
    binlog(convert(UInt32, v))

efficient log_2 of integers
https://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
"""
binlog(v::Integer) = binlog(convert(UInt32, v))
function binlog(v::UInt32)
    r =     (v > 0xFFFF) << 4; v >>= r;
    shift = (v > 0xFF  ) << 3; v >>= shift; r |= shift;
    shift = (v > 0xF   ) << 2; v >>= shift; r |= shift;
    shift = (v > 0x3   ) << 1; v >>= shift; r |= shift;
                                            r |= (v >> 1);
end

"""
    Faber(n)

change of indexation: from single indexation n to double indexation i,j
"""
function Faber(n)
    return (binlog(n)), Int(n - 1<<(binlog(n)))
end

"""
    OrderBasis(i, j)

compare two index (single indexation) and order it from the smaller
"""
function OrderBasis(i, j)
    if i <= j
        return (i, j)
    else
        return (j, i)
    end
end

"""
    OrderBasis(i1, j1, i2, j2)

compare two index (double indexation) and order it from the smaller
"""
function OrderBasis(i1 ,j1 ,i2 ,j2)
    if i1 <= i2
        return (i1, j1, i2, j2)
    else
        return (i2, j2, i1, j1)
    end
end
"""
    ShareDomain(i,j)

check if the basis relative to index i and index j (single indexation) share the same support
return 1 if they do
return 0 if they don't
"""
function ShareDomain(i,j)
    if i == j return 1
    else
        (i, j) = OrderBasis(i, j)
        i1, j1 = Faber(i)
        i2, j2 = Faber(j)
        if j1*2<<(i2 - i1 - 1) <= j2 < (j1 + 1)*2<<(i2 - i1-1)
            return 1
        else
            return 0
        end
    end
end

"""
    SharedSupport(n, L)

Collect all the indexes (single indexation) relative to the basis function sharing
sharing the same support of the function ϕ_n up to level L
"""
function SharedSupport(n,L)
    a = Int64[]
    for i in 1:(2<<(L)) - 1
        if ShareDomain(n,i) == 1 # && n != i
            push!(a, i)
        end
    end
    return a
end

"""
    SharedSupport(i,j, L)

Collect all the indexes (with single indexation) relative to the basis functions
sharing the same support of the function ϕ_{i,j} up to  level L
"""
function SharedSupport(i, j, L)
    n = Faber(i,j)
    a = SharedSupport(n, L)
    return a
end
"""
    NotSharedSupport(i,j, L)

Collect all the indexes (with single indexation) relative to the basis functions
which do the same support of the function ϕ_{i,j} up to  level L
"""
function NotSharedSupport(n, L)
    a = Int64[]
    for i in 1:(2<<(L)) - 1
        if ShareDomain(n,i) == 0
            push!(a, i)
        end
    end
    return a
end

"""
    FS

the strcut
```
struct Fs
    i::Int16    #idenx i (double indexation)
    j::Int16    #index j (double indexation)
    n::Int64    #idenx n (single indexation)
    lb::Float64 #pper bound domain
    ub::Float64 #lower Bound domain [not used?]
    range::Float64 #ub - lb
    nhb::Array{Int64} #neightbourhoods
    nhb::Array{Int64} #not neightbourhoods
    supr::Float64 #suprimum(ϕ)
    δ::Float64
end
```
"""
struct Fs
    i::Int16    #idenx i (double indexation)
    j::Int16    #index j (double indexation)
    n::Int64    #idenx n (single indexation)
    lb::Float64 #Upper bound domain
    ub::Float64 #Lower Bound domain [not used?]
    range::Float64
    nhb::Array{Int64} #neightbourhoods
    notnhb::Array{Int64}
    supr::Float64 #suprimum(ϕ)
    δ::Float64 # range*suprimum(ϕ)
    function Fs(i, j, L, T) #inizialization
        n = Faber(i, j)
        new(i, j, n, T*j/2^i, T*(j+1)/2^i, T/2^i, SharedSupport(n,L), NotSharedSupport(n,L), 2^(-i/2 - 1)*sqrt(T), T/2^i*2^(-i/2 -1)*sqrt(T))
    end
    function Fs(n, L, T) #inizialization
        i, j = Faber(n)
        Fs(i, j, L, T)
    end
end





"""
    generate(L, T)

generate all the set of basis used consistenly
"""
function generate(L, T)
    ϕ = Array{Fs}(undef, 2<<L - 1)
    for  n in 1:2<<L -1
        ϕ[n] = Fs(n, L, T)
    end
    ϕ
end

"""

"""
function Λ(x, T)
    return max(0.0, sqrt(T)/2 - abs(x/sqrt(T) - sqrt(T)/2))
end
#T = 3
#plot([Λ(x, T) for x=0:0.001:T])
"""


"""
function Λ(x, l, k, T)
    x = mod(x, T)   #ok
    2^(-l/2)*Λ((x)*(1<<l) - k*T, T)
end
