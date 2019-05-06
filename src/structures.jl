#ordering from the larger basis to the smaller
function OrderBasis(i, j)
    if i <= j
        return (i, j)
    else
        return (j, i)
    end
end

function OrderBasis(i1 ,j1 ,i2 ,j2)
    if i1 <= i2
        return (i1, j1, i2, j2)
    else
        return (i2, j2, i1, j1)
    end
end

function ShareDomain(i,j)
    if i == j return 1
    else
        (i, j) = OrderBasis(i, j)
        i1, j1 = Faber(i)
        i2, j2 = Faber(j)
        if j1*2^(i2 - i1) <= j2 < (j1 + 1)*2^(i2 - i1)
            return 1
        else
            return 0
        end
    end
end
function SharedSupport(n,L)
    a = Int64[]
    for i= 1:(1<<(L+1)) - 1
        if ShareDomain(n,i) == 1 && n != i
            push!(a, i)
        end
    end
    return a
end

function SharedSupport(i,j,L)
    n = Faber(i,j)
    a = SharedSupport(n, L)
    return a
end

struct GFSbase
    i::Int16
    j::Int16    #identifiers 2 i
    n::Int64    #identifier 1 i
    L::Int64    #max level
    lb::Float64 #Upper bound domain
    ub::Float64 #Lower Bound domain
    range::Float64
    nhb::Array{Int64} #neightbourhood
    T::Float64
    δ::Float64 # range*suprimum(ϕ)
    function GFSbase(i, j, L, T) #inizialization
        n = Faber(i, j)
        new(i, j, n, L , T*j/2^i, T*(j+1)/2^i, T/2^i, SharedSupport(n,L), T, T/2^i*2^(-i/2 -1)*sqrt(T))
    end
    function GFSbase(n, L, T) #inizialization
        i, j = Faber(n)
        GFSbase(i, j, L, T)
    end
end

struct Skeleton
    t::Float64
    ξ::Array{Float64, 1}
end
