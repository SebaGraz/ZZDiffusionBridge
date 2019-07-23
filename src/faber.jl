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
    Λ(x, T)

function needed for construction of any Faber Schauder function i,j
see   Λ(x, l, k, T). This function should not be visible.
"""
function Λ(x, T)
    return max(0.0, sqrt(T)/2 - abs(x/sqrt(T) - sqrt(T)/2))
end

"""
    Λ(x, l, k, T)

Faber Schauder function `l`, `k` (double indexation) construct as in definition
rescaled by length `T` and evaluated at `x`
"""
function Λ(x, l, k, T)
    x = mod(x, T)   #ok
    2^(-l/2)*Λ((x)*(1<<l) - k*T, T)
end

"""
    IntQuad1(a, b, n1, n2)
integral of quadratic function1 with `a` quadratic term
`b` linear term and`n1` `n2` uppper and lower bound:

∫_n2^n1 a*t^2 + b*t dt
"""
function IntQuad1(a, b, n1, n2)
    a*(n1^3 - n2^3)/3 + b*(n1^2 - n2^2)/2
end


"""
    flip(i,j)

flip index j to the left to ease the integration:
int_0^1 psi_i (t) psi_j(t) dt

ex flip(1,1) = 0 flip(1,0) = 0 (do nothing), flip(3,5)= 2
"""
function flip(i,j)
    bound = ((1/2)*2^i - 1/2)
    if bound < j
        j = 2*bound - j
    end
    return Integer(j)
end
"""
    flip(i,j)

flip index j on the other side to ease the computation
of bound2(n,T)

ex. flip2(1,1) = 0, flip2(1,0) = 1, flip2(2,0) = 3, flip2(2,1) = 2.

"""
function flip2(i,j)
    bound = ((1/2)*2^i - 1/2)
    j = 2*bound - j
    return Integer(j)
end


"""
    int_prod_fs(i,j, T)

compute Φ_{ij} = ∫_0^T ϕ_i ϕ_j dt
ϕ_i being the ith (single indexation) Faber Schauder function
"""
#int_0^1 psi_i (t) psi_j(t) dt with i<j
function int_prod_fs(i,j, T)
    (i, j) = OrderBasis(i, j)
    i1, j1 = Faber(i) #ordering  i < j --> i1 < i2
    i2, j2 = Faber(j)
    i3 = i2 - i1
    j3 = j2 - 2^(i2-i1)*j1
    j3 = flip(i3,j3)
    A = j3/2^i3
    B = (j3 + 1/2)/2^i3
    C = (j3 + 1)/ 2^i3
    return T^2*2^((i2 - 5*i1)/2)*(IntQuad1(1,-j3/2^i3, B, A) - IntQuad1(1,-(j3+1)/2^i3, C, B))
end

"""
    intfs(i, T)

compute Φ_i = ∫_0^T ϕ_{ij} ds
"""
function intfs(n, T)
    i, j = Faber(n)
    return T^(1.5)*2^(-(3i + 4)/2 )
end

"""
    generate_matrix(L, T)
generate vector whose element is ∫_0^T ϕ(i1,j1)*ϕ(i2,j2) dt
for each i1,j1,i2,j2 until truncation level `L`
"""

function generate_matrix(L::Int64, T)
    dim = 2^(L+1) -1
    V = zeros(dim,dim)
    for i in 1:dim
        for j in 1:i
            if i == j
                i1, j1 = Faber(i)
                V[i,i] = T^2* 2.0^(-2*i1 - 2)/3
            else
                if ShareDomain(i,j)==1
                    res = int_prod_fs(i, j, T)
                    V[i,j], V[j,i] = (res, res)
                else
                    V[i,j], V[j,i] = (0 , 0)
                end
            end
        end
    end
    return V
end


"""
    generate_vector(L, T)
generate vector whose element is ∫_0^T ϕ(i,j) dt
for each i,j until truncation level `L`
"""

function generate_vector(L, T)
    [intfs(i, T) for i in 1:(2<<L -1)]
end

"""
    bound2(L, T)
computes ∫_0^T ϕ(n)(t/T) dt
N.B. `n` single inexation
"""
function bound1(n, T)
    i,j = Faber(n)
    A = j / (2^i)
    B = (0.5 + j)/ (2^i)
    C = (1.0 + j)/ (2^i)
    return T^(3/2)*2^(-i/2)*(IntQuad1(2^i,-j, B, A) + IntQuad1(-2^i, j + 1, C, B))
end

"""
    bound2(L, T)
computes ∫_0^T ϕ(n)(1 - t/T) dt
N.B. `n` single inexation
"""
function bound2(n, T)
    i,j = Faber(n)
    j = flip2(i,j)
    A = j / (2^i)
    B = (0.5 + j)/ (2^i)
    C = (1.0 + j)/ (2^i)
    return T^(3/2)*2^(-i/2)*(IntQuad1(2^i,-j, B, A) + IntQuad1(-2^i, j + 1, C, B))
end

"""
    generate_bound1(L, T)
generate vector whose element is ∫_0^T ϕ(i,j)(t/T) dt
for each i,j until truncation level `L`
"""
function generate_bound1(L, T)
    [bound1(n,T) for n in 1:1:2<<L-1]
end

"""
    generate_bound2(L, T)
generate vector whose element is ∫_0^T ϕ(i,j)(1 - t/T) dt
for each i,j until truncation level `L`
"""
function generate_bound2(L, T)
    [bound2(n,T) for n in 1:1:2<<L-1]
end

"""
    MCintegration(k::Fs)

take a random point inside the support of the basis function k
"""
function MCintegration(k::Fs)
    rand()*k.range + k.lb
end

"""
    Sol2E(a::Float64, b::Float64, c::Float64)

Solve a quadratic equation.
TODO add errors when the solution is not real.
solve quadratic equation ax^2 + bx + c = 0
"""
function Sol2E(a::Float64, b::Float64, c::Float64)
    return (-b + sqrt(b^2 - 4a*c))/2a
end
