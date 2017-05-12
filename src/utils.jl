function relorientation(face, simplex)

    v = setdiff(simplex, face)
    length(v) == 1 || return 0

    # find the position of the missing vertex
    v = v[1]
    i = Base.findfirst(simplex, v)
    s = (-1)^(i-1)

    # remove that vertex from the simplex
    face2 = Array{Int}(length(simplex)-1)
    for j in 1 : i-1
        face2[j] = simplex[j]
    end
    for j in i : length(simplex)-1
        face2[j] = simplex[j+1]
    end

    # get the permutation that maps face to face2
    #p = indexin(face, face2)
    p = [ findfirst(face2,v) for v in face ]

    return s * levicivita(p) * i
end
