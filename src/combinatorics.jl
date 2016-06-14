export getcommonedge


const levicivita_lut = cat(3,
    [0 0  0;  0 0 1; 0 -1 0],
    [0 0 -1;  0 0 0; 1  0 0],
    [0 1  0; -1 0 0; 0  0 0])



# Levi-Civita symbol of a permutation.
# The parity is computed by using the fact that a permutation is odd if and
# only if the number of even-length cycles is odd.
# Returns 1 is the permutarion is even, -1 if it is odd and 0 otherwise.
function levicivita(p)
    n = length(p)

    if n == 3
        @inbounds valid = (0 < p[1] <= 3) * (0 < p[2] <= 3) * (0 < p[3] <= 3)
        return valid ? levicivita_lut[p[1], p[2], p[3]] : 0
    end

    todo = trues(n)
    first = 1
    cycles = flips = 0

    while cycles + flips < n
        first = findnext(todo, first)
        (todo[first] $= true) && return 0
        j = p[first]
        (0 < j <= n) || return 0
        cycles += 1
        while j ≠ first
            (todo[j] $= true) && return 0
            j = p[j]
            (0 < j <= n) || return 0
            flips += 1
        end
    end

    return iseven(flips) ? 1 : -1
end



import Base.indexin
@inline indexin(V::FixedArray, A::AbstractArray) = [ findfirst(A,v) for v in V ]



# given a simplex and a face returns:
# +v if face is the v-th face of the simplex oriented according to the simplex
# -v if face is the v-th face of the simplex oriented oppositely to the simplex
# 0 is face is not a face of the simplex
function relorientation(face, simplex)

    v = setdiff(simplex, face)
    # if length(v) != 1
    #     return 0
    # end
    length(v) == 1 || return 0

    # find the position of the missing vertex
    v = v[1]
    i = Base.findfirst(simplex, v)
    s = (-1)^(i-1)

    # remove that vertex from the simplex
    face2 = Array(Int, length(simplex)-1)
    for j in 1 : i-1
        face2[j] = simplex[j]
    end
    for j in i : length(simplex)-1
        face2[j] = simplex[j+1]
    end

    # get the permutation that maps face to face2
    #p = findfirst(face2,face)
    #p = findfirst(face2, face)
    p = indexin(face, face2)

    return s * levicivita(p) * i
end


"""
    getcommonedge(cell1, cell2) -> e1, e2, edge

Returns in edge the common vertices of cell1 and cell2. e1 contains the index
of the vertex of cell1 opposite to this common edge, and with a plus or minus
sign depending on whether the orientation of the common edge is along or
against the internal orientation of cell1. Similar for e2.
"""
function getcommonedge(cell1, cell2)
    isct = intersect(cell1, cell2)
    relorientation(isct, cell1), relorientation(isct, cell2), isct
end
