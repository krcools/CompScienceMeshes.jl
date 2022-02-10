"""
gives the ordering for nedelec2 (divergence)
"""
function relorientation_gen(face::SArray{Tuple{3},T,1,3}, tet::SArray{Tuple{4},T,1,4}) where {T}
    # v = setdiff(tet,face)
    # length(v) == 1 || return 0

    # a = something(findfirst(isequal(v[1]), tet),0)

    a = findfirst(x -> !(x in face), tet)
    a == nothing && return 0

    @assert 1 <= a <= 4
    q = findnext(x ->!(x in face), tet, a+1)
    q == nothing || return 0

    v = tet[a]

    Face = deleteat(tet,a)
    # @show a tet face Face indexin(face,Face)
    p = parity(Vector{Int}(indexin(face,Face)))
    return a * (-1)^a * (-1)^p

    # face = [v[1],face[1],face[2],face[3]]
    # w = sortperm(face)
    # b = parity(w)
    # return -a*(-1)^b
end

const RelOp34 = cat(
    [0 0 0 0; 0 0 4 -3; 0 -4 0 2; 0 3 -2 0],
    [0 0 -4 3; 0 0 0 0; 4 0 0 -1; -3 0 1 0],
    [0 4 0 -2; -4 0 0 1; 0 0 0 0; 2 -1 0 0],
    [0 -3 2 0; 3 0 -1 0; -2 1 0 0; 0 0 0 0],
    dims=3)

function relorientation(face::SArray{Tuple{3},T,1,3}, tet::SArray{Tuple{4},T,1,4}) where {T}

    i = findfirst(==(face[1]),tet); i == nothing && return 0
    j = findfirst(==(face[2]),tet); j == nothing && return 0
    k = findfirst(==(face[3]),tet); k == nothing && return 0
    return RelOp34[i,j,k]

end

"""
gives the ordering for nedelec1 (curl)
if a tetrahedron is given by simplex(a,b,c,d), the order is:
    a-b
    a-c
    a-d
    b-c
    b-d
    c-d
"""

const _relorient24 = @SMatrix[
    0 1 2 3
    -1 0 4 5
    -2 -4 0 6
    -3 -5 -6 0]

# function relorientation(edge::SArray{Tuple{2},Int64,1,2}, tet::SArray{Tuple{4},Int64,1,4})
@inline function relorientation(edge::SVector{2,Int}, tet::SVector{4,Int})

    # relorient24 = [
    # 0 1 2 3
    # -1 0 4 5
    # -2 -4 0 6
    # -3 -5 -6 0]

    # # println("relorientation 2-4")
    # e1 = edge[1]
    # a = 0
    # for v in tet
    #     a = a + 1
    #     if e1 == v
    #         break
    #     end
    # end
    # a == 5 && return 0

    # e2 = edge[2]
    # b = 0
    # for v in tet
    #     b = b + 1
    #     if e2 == v
    #         break
    #     end
    # end
    # b == 5 && return 0
    # return relorient24[a,b]

    a = findfirst(==(edge[1]), tet)
    a == nothing && return 0

    b = findfirst(==(edge[2]), tet)
    b == nothing && return 0

    return _relorient24[a,b]

end

# function relorientation(edge::SArray{Tuple{2},Int64,1,2}, tet::SArray{Tuple{4},Int64,1,4})
#     # v = setdiff(tet, edge)
#     # length(v) == 2 || return 0

#     count(x -> !(x in edge), tet) != 2 && return 0

#     w1 = findfirst(isequal(edge[1]), tet)
#     w2 = findfirst(isequal(edge[2]), tet)

#     @assert w1 != nothing
#     @assert w2 != nothing

#     t = tetschoice(w1,w2)
#     s = sign(w2-w1)

#     return s*t
# end

const RelOp23 = [0 3 -2; -3 0 1; 2 -1 0]
function relorientation(face::SVector{2,Int}, cell::SVector{3,Int})
    i = findfirst(==(face[1]), cell); i == nothing && return 0
    j = findfirst(==(face[2]), cell); j == nothing && return 0
    return RelOp23[i,j]
end


function tetschoice(w1::Int64,w2::Int64)
    a = min(w1,w2)
    b = max(w1,w2)
    if a == 1
        return b-1
    elseif a == 2
        return b+1
    elseif a == 3
        return b+2
    end
end

function relorientation(node::SArray{Tuple{1},Int64,1,1}, face::SArray{Tuple{N},Int64,1,N} where {N})
    idx = something(findfirst(isequal(node[1]), face), 0)
    sgn = (-1)^(idx-1)
    return sgn * idx
end

function relorientation(face, simplex)

    v = setdiff(simplex, face)
    length(v) == 1 || return 0

    # find the position of the missing vertex
    v = v[1]
    i = something(findfirst(isequal(v), simplex),0)
    s = (-1)^(i-1)

    # remove that vertex from the simplex
    face2 = Array{Int}(undef,length(simplex)-1)
    for j in 1 : i-1
        face2[j] = simplex[j]
    end
    for j in i : length(simplex)-1
        face2[j] = simplex[j+1]
    end

    # get the permutation that maps face to face2
    p = [something(findfirst(isequal(v),face2),0) for v in face]

    return s * levicivita(p) * i
end

relorientation(p::SVector{4}, q::SVector{3}) = relorientation(q,p)
relorientation(p::SVector{4}, q::SVector{2}) = relorientation(q,p)
relorientation(p::SVector{4}, q::SVector{1}) = relorientation(q,p)
relorientation(p::SVector{3}, q::SVector{2}) = relorientation(q,p)
relorientation(p::SVector{3}, q::SVector{1}) = relorientation(q,p)
relorientation(p::SVector{2}, q::SVector{1}) = relorientation(q,p)