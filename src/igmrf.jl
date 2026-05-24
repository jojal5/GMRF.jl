struct iGMRF
    G::GridStructure
    rankDeficiency::Int64
    κ::Float64              # Precision of the field
end

function Base.show(io::IO, obj::iGMRF)

    println(io, "iGMRF")
    println(io, "G :")
    showGridStructure(io, obj.G, prefix = "\t\t\t")
    println(io)
    println(io, "rankDeficiency :\t", obj.rankDeficiency)
    println(io, "κ :\t\t\t", obj.κ)

end

function iGMRF(m₁::Integer, m₂::Integer, order::Integer, κ::Real)::iGMRF

    # Gives the adjacency matrix W for the iGMRF of order 1 or 2 on the regular
    # grid of size (m1 * m2).

    @assert order == 1 || order == 2 "the order should be either 1 or 2."

    if order == 1

        nbs, W = first_order_lattice_neighbors(m₁, m₂)
        condIndSubset = fo_condindsubsets(m₁, m₂)
        rankdef = 1

    else

        nbs, W = so_nbs(m₁, m₂)
        condIndSubset = so_condindsubsets(m₁, m₂)
        rankdef = 3

    end

    W̄ = W - spdiagm(length.(nbs))

    G = GridStructure((m₁, m₂), nbs, condIndSubset, W, W̄)

    return iGMRF(G, rankdef, κ)

end




function fo_condindsubsets(m₁::Integer, m₂::Integer)::Vector{Vector{Integer}}


    condIndSubsetIndex = 2*ones(Int64,m₁,m₂)
    condIndSubsetIndex[1:2:end,1:2:end] .= 1
    condIndSubsetIndex[2:2:end,2:2:end] .= 1

    return Array[findall(vec(condIndSubsetIndex) .==i) for i=1:2]

end

function so_condindsubsets(m₁::Integer, m₂::Integer)::Vector{Vector{Integer}}

    condIndSubsetIndex = zeros(Int64,m₁,m₂)

    condIndSubsetIndex[1:3:end,1:4:end] .= 1
    condIndSubsetIndex[2:3:end,3:4:end] .= 1

    condIndSubsetIndex[1:3:end,2:4:end] .= 2
    condIndSubsetIndex[2:3:end,4:4:end] .= 2

    condIndSubsetIndex[1:3:end,3:4:end] .= 3
    condIndSubsetIndex[2:3:end,1:4:end] .= 3

    condIndSubsetIndex[1:3:end,4:4:end] .= 4
    condIndSubsetIndex[2:3:end,2:4:end] .= 4

    condIndSubsetIndex[3:3:end,1:3:end] .= 5
    condIndSubsetIndex[3:3:end,2:3:end] .= 6
    condIndSubsetIndex[3:3:end,3:3:end] .= 7

    return Array[findall(vec(condIndSubsetIndex) .==i) for i=1:7]

end

function rand(F::iGMRF)::Vector{<:Real}

    @assert F.rankDeficiency == 1 || F.rankDeficiency == 3 "The rank deficiency should be either 1 or 3"

    κ = F.κ
    W = F.G.W
    m₁ = F.G.gridSize[1]
    m₂ = F.G.gridSize[2]
    m = m₁ * m₂

    if F.rankDeficiency == 1

        e₁ = ones(m,1)

        A = e₁

        Q = κ*W + e₁*e₁'

    else

        e₁ = ones(m)
        e₂ = repeat(1:m₁, m₂)
        e₃ = repeat(1:m₂,inner = m₁)

        A = hcat(e₁,e₂,e₃)

        Q = κ*W + e₁*e₁' + e₂*e₂' + e₃*e₃'

    end

    C = cholesky(Q)
    L = C.L

    z = randn(m)

    x = L'\z

#     V = zeros(m,size(A,2))
#     for ii=1:size(A,2)
#        V[:,ii] = C\A[:,ii]
#     end
    V = C\A
    W = A'*V
    U = W\(V')
    c = A' * x
    y = x - U' * c

    return y

end

function logpdf(F::iGMRF, y::Array{<:Real})::Real

    κ = F.κ

    W = F.G.W
    m = F.G.gridSize[1] * F.G.gridSize[2]

    k = F.rankDeficiency

    v = κ*W*y
    q = y'*v

    lpdf =  .5*(m-k)*log(κ) - .5*q

    return lpdf

end

function fullconditionals(F::iGMRF, y::Vector{<:Real})::Vector{NormalCanon}

    κ = F.κ

    W̄ = F.G.W̄
    W = F.G.W

    Q = κ * Array(diag(F.G.W))
    h = -κ*(W̄*y)

    pd = NormalCanon.(h,Q)

    return pd

end

function fullcondlogpdf(F::iGMRF, y::Vector{<:Real})::Vector{<:Real}

    pd = fullconditionals(F::iGMRF,y::Vector{<:Real})

    clpdf = logpdf.(pd,y)

    return clpdf

end

function getconditional(F::GMRF.iGMRF, B::Vector{<:Integer}, x::Vector{<:Real})::MvNormalCanon

    W = F.G.W

    sort!(B)

    A = setdiff(1:(F.G.gridSize[1] * F.G.gridSize[2]), B)

    Waa = W[A,A]
    Wab = W[A,B]

    h = -Wab*x*F.κ

    J = Array(F.κ*Waa)

    pd = MvNormalCanon(h,J)

    return pd

end
