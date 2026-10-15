# N: the number of nodes (aka the dimension plus one)
# P: the order of the simplex
# Q: the number of vertices (depends on N and the order of the simplex)
struct HOSimplexGraph{N,P,Q} <: AbstractGraph
    indices::SVector{Q,Int}
end

function HOSimplexGraph{N,P}(idcs::Vararg{Int,Q}) where {N,P,Q}
    @assert Q == binomial(N+P-1,P)
    HOSimplexGraph{N,P,Q}(SVector{Q,Int}(idcs))
end

# edges of a triangle
function faces(g::HOSimplexGraph{3,P,Q}, ::Type{Val{1}}) where {P,Q}

    N=3

    Nout = N-1
    Pout = P

    F = HOSimplexGraph{Nout,Pout}
    I = g.indices


    (
        F(I[2], I[3+P:3+2P-2]..., I[3]),
        F(I[3], I[3+2P-1:3+3P-3]..., I[1]),
        F(I[1], I[3+1:3+P-1]..., I[2]),
    )
end


