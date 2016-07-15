
#import Base.getindex

export simplex
export dimension, universedimension, vertextype
export vertices, tangents, volume
export barytocart, carttobary, centroid

export FlatCellNM


# U: the dimension of the universe
# D: the dimension of the manifold
# N: the number of vertices
# T: the type of the coordinates
# C: the complimentary dimension (should always be U-D)
immutable FlatCellNM{U,D,C,N,T}
    vertices::Vec{N,Vec{U,T}}
    tangents::Vec{D,Vec{U,T}}
    normals::Vec{C,Vec{U,T}}
    volume::T
end




"""
    coordtype(simplex)

Return coordinate type used by simplex.
"""
coordtype{U,D,C,N,T}(::Type{FlatCellNM{U,D,C,N,T}}) = T
coordtype(p::FlatCellNM) = coordtype(typeof(p))

"""
    vertices(simplex)

Return an indexable iterable to the vertices defining simplex
"""
vertices(p::FlatCellNM) = p.vertices



"""
    tangents(simplex, i)

Return the i-th tangent to simplex.
"""
tangents(p::FlatCellNM,i) = p.tangents[i]



"""
    normal(simplex)

Return the unit normal to a hypersimplex (i.e. dim(simplex) =
dim(universe)-1). Undefined for non-hyper simplices.
"""
normal(p::FlatCellNM) = p.normals[1]

"""
    volume(simplex)

Return the volume of the simplex.
"""
volume(p::FlatCellNM) = p.volume


"""
    dimension(simplex)

Return the manifold dimension of the simplex.
"""
#dimension(p::FlatCellNM) = length(p.vertices) - 1
dimension{U,D,C,N,T}(::Type{FlatCellNM{U,D,C,N,T}}) = D
dimension(p::FlatCellNM) = dimension(typeof(p))


"""
  universedimension(p)

Return the dimension of the universe in which `p` is embedded.
"""
universedimension{U,D,C,N,T}(::Type{FlatCellNM{U,D,C,N,T}}) = U
universedimension(p::FlatCellNM) = universedimension(typeof(p))


"""
    getindex(simplex, I)
    simplex[I]

Get the vertices at index I (scalar or array) defining the simplex
"""
getindex(p::FlatCellNM, I::Union{Number,Vec,Array}) = p.vertices[I]



"""
    simplex(vertices)
    simplex(v1, v2, ...)
    simplex(vertices, Val{D})

Build a D-dimensional simplex. The vertices can be passed in
an array (static or dynamic), or supplied separately. If the
length of the array is not part of its type, the speed of the
construction can be improved by supplying an extra Val{D}
argument. In case it is not clear from the context whether
the vertex array is dynamically or statically sized, use the
third form as it will not incur notable performance hits.

Note that D is the dimension of the simplex, i.e. the number
of vertices supplied minus one.
"""
@generated function simplex{D1,P}(vertices::Vec{D1,P})
    U = length(P)
    D = D1 - 1
    C = U-D
    T = eltype(P)
    xp1 =:(())
    for i in 1:D
        push!(xp1.args, :(vertices[$i]-vertices[end]))
    end
    xp2 = :(Vec{$D,P}($xp1))
    quote
        tangents = $xp2
        normals, volume = _normals(tangents, Val{$C})
        FlatCellNM(vertices, tangents, normals, volume)
    end
    # tangents = Vec{D1-1,P}((v[1]-v[end], v[2]-v[end], ...))
    # normals, volume = _normals(tangents, Val{length(P)-D1+1})
    # FlatCellNM(vertices, tangents, normals, volume)
end

simplex(vertices...) = simplex(Vec((vertices...)))

@generated function simplex{D}(vertices, ::Type{Val{D}})
    P = eltype(vertices)
    D1 = D + 1
    xp = :(())
    for i in 1:D1
        push!(xp.args, :(vertices[$i]))
    end
    :(simplex(Vec{$D1,$P}($xp)))
end


function _normals(tangents, ::Type{Val{1}})
    PT = eltype(tangents)
    D  = length(tangents)
    T  = eltype(PT)

    n = zeros(T,D+1)
    b = Array(T,D,D)

    for i in 1:D+1
        fill!(b, zero(T))
        for j in 1:D
            for k in 1:i-1
                b[k,j] = tangents[j][k]
            end
            for k in i:D
                b[k,j] = tangents[j][k+1]
            end
        end
        n[i] = (-1)^(i-1) * det(b)
    end

    n *= (-1)^D / norm(n)
    normals = Vec{1,PT}([PT(n)])

    metric = T[dot(tangents[i], tangents[j]) for i in 1:D, j in 1:D]
    volume = sqrt(abs(det(metric))) /  D

    return normals, volume

end




function _normals{C}(tangents, ::Type{Val{C}})
    PT = eltype(tangents)
    D  = length(tangents)
    U = length(PT)
    T  = eltype(PT)

    metric = T[dot(tangents[i], tangents[j]) for i in 1:D, j in 1:D]
    volume = sqrt(abs(det(metric))) / D

    # Fix this. This function needs to become gneerated
    normals = Vec{C,PT}([zero(PT) for i in 1:C])

    return normals, volume
end



"""
    barytocart(simplex, uv)

Returns the point in the simplex with barycentric coordinates uv
"""
function barytocart(mani::FlatCellNM, u)
    r = last(mani.vertices)
    for i in 1 : dimension(mani)
        ti = mani.tangents[i]
        ui = u[i]
        #r += mani.tangents[i] * u[i]
        r += ti * ui
    end
    return r
end



function carttobary{U,D,C,N,T}(p::FlatCellNM{U,D,C,N,T}, cart)

    G = [dot(p.tangents[i], p.tangents[j]) for i in 1:D, j in 1:D]

    w = [dot(p.tangents[i], cart - p.vertices[end]) for i in 1:D]
    u = G \ w

    return Vec(u)
end
