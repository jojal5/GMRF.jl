
function structure_igmrf(m₁::Integer,m₂::Integer,order::Integer,κ::Real)

#=Gives the adjacency matrix W for the iGMRF of order 1 or 2 on the regular
grid of size (m1 * m2). =#

    m = m₁*m₂

    if order == 1

        rankDeficiency = 1

        # 1-off diagonal elements
        v = ones(Int64,m₁)
        v[end] = 0
        V = repeat(v,outer=m₂)
        pop!(V)

        # n-off diagonal elements
        U = ones(Int64,m₁*(m₂-1))

        # get the upper triangular part of the matrix
        D = sparse(1:(m-1), 2:m, V, m, m) + sparse(1:(m-m₁),(m₁+1):m, U, m, m)

        # make W symmetric
        D = Symmetric(D,:U)

        # Compute the list of neighbors for each node
        nbs =  Array{Int64,1}[]
        for i = 1:m
            push!(nbs,findall(D[:,i] .> 0))
        end

        nnbs = length.(nbs)

        # Put the number of neighbors on the diagonal
        W = -D + spdiagm(0 => nnbs)


    elseif order==2

        rankDeficiency = 3

        # Alternative by adding molecules. There should not be missing values in the grid.

        W = spzeros(Int64,m,m)
        pos = reshape(1:m,m₁,m₂)

        for i=1:m₁
            for j=1:m₂

                S = zeros(Int64,m₁,m₂)

                if (i-2>0)
                   S[i-2:i,j] =  S[i-2:i,j] + [1, -2, 1]
                end

                if (i+2<=m₁)
                   S[i:i+2,j] =  S[i:i+2,j] + [1, -2, 1]
                end

                if (j-2>0)
                   S[i,j-2:j] =  S[i,j-2:j] + [1,-2, 1]
                end

                if (j+2<=m₂)
                   S[i,j:j+2] =  S[i,j:j+2] + [1,-2, 1]
                end



                if (i-1>0) && (i+1<=m₁)
                    S[i-1:i+1,j] = S[i-1:i+1,j] + [-2, 4, -2]
                end

                if (j-1>0) && (j+1<=m₂)
                    S[i,j-1:j+1] = S[i,j-1:j+1] + [-2, 4, -2]
                end



                if (i-1>0) && (j+1<=m₂)
                    S[i-1:i,j:j+1] = S[i-1:i,j:j+1] + [-2 2; 2 -2]
                end

                if (i+1<=m₁) && (j+1<=m₂)
                    S[i:i+1,j:j+1] = S[i:i+1,j:j+1] + [2 -2; -2 2]
                end

                if (i-1>0) && (j-1>0)
                    S[i-1:i,j-1:j] = S[i-1:i,j-1:j] + [2 -2; -2 2]
                end

                if (i+1<=m₁) && (j-1>0)
                    S[i:i+1,j-1:j] = S[i:i+1,j-1:j] + [-2 2; 2 -2]
                end

                W[:,pos[i,j]] = S[:]

            end
        end

        # Compute the list of neighbors for each node
        nbs =  Array{Int64,1}[]
        for i = 1:m
            push!(nbs,findall(W[:,i] .< 0))
        end

        nnbs = length.(nbs)

    end

    condIndSubset = get_condindsubsets(m₁,m₂,order)


    G = GridStructure(m, (m₁,m₂), nbs, nnbs, condIndSubset)

    W̄ = W - sparse(diagm(nnbs))

    F = iGMRF(G, rankDeficiency, κ, W, W̄)

    return F

end

function get_condindsubsets(m₁::Integer,m₂::Integer,order::Integer)

    if order == 1

        condIndSubsetIndex = 2*ones(Int64,m₁,m₂)
        condIndSubsetIndex[1:2:end,1:2:end] .= 1
        condIndSubsetIndex[2:2:end,2:2:end] .= 1

        condIndSubset = Array[findall(vec(condIndSubsetIndex) .==i) for i=1:2]

    elseif order == 2

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

        condIndSubset = Array[findall(vec(condIndSubsetIndex) .==i) for i=1:7]

    end

    return condIndSubset

end

function rand(F::iGMRF)

    κ = F.κ
    W = F.W
    m₁ = F.G.gridSize[1]
    m₂ = F.G.gridSize[2]
    m = F.G.m

    if F.rankDeficiency == 1

        e₁ = ones(m,1)

        A = e₁

        Q = κ*W + e₁*e₁'

    elseif F.rankDeficiency == 3

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

function logpdf(F::iGMRF,y::Array{Float64})

    κ = F.κ

    W = F.W
    m = F.G.m

    k = F.rankDeficiency

    v = κ*W*y
    q = y'*v

    lpdf =  .5*(m-k)*log(κ) - .5*q

    return lpdf

end

function fullconditionals(F::iGMRF,y::Vector{<:Real})

    κ = F.κ

    W̄ = F.W̄
    W = F.W

    Q = κ * Array(diag(F.W))
    h = -κ*(W̄*y)

    pd = NormalCanon.(h,Q)

    return pd

end

function fullcondlogpdf(F::iGMRF,y::Vector{<:Real})

    pd = fullconditionals(F::iGMRF,y::Vector{<:Real})

    clpdf = logpdf.(pd,y)

    return clpdf

end
