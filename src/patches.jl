
import Base.getindex

# Factory methods
export simplex

# Interface for any Patch
export barytocart, vertices, dimension, volume, carttobary

# Interface for FlatPatch
export tangents, dimension, centroid
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
coordtype(p::FlatCellNM) = eltype(eltype(p.vertices))

"""
    vertices(simplex)

Return an indexable iterable to the vertices defining simplex
"""
vertices(p::FlatCellNM) = p.vertices

"""
    tangents(simplex)

Return an indexable iterable to the tangent vectors to simplex
"""
tangents(p::FlatCellNM) = p.tangents

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
dimension(p::FlatCellNM) = length(p.vertices) - 1

"""
    getindex(simplex, I)
    simplex[I]

Get the vertices at index I (scalar or array) defining the simplex
"""
getindex(p::FlatCellNM, I) = p.vertices[I]



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

# This one is, as expected quote a bit slower...
simplex(vertices::Array) = simplex(Vec(vertices))


# function patch{P<:FixedArray,D}(verts::Vector{P}, ::Type{Val{D}})
#   @assert length(verts) == D+1
#   U = length(P)
#   T = eltype(P)
#   # TODO: to this in a generated function
#   vertices = Vec{D+1,P}(verts)
#   tangents = Vec{D,P}([verts[i]-verts[end] for i in 1:D])
#   normals, volume = _normals(tangents, Val{U-D})
#   FlatCellNM(vertices, tangents, normals, volume)
# end
#
# function patch{T}(verts::Vector{Vec{3,T}}, ::Type{Val{2}})
#   @assert length(verts) == 3
#   PT = eltype(verts)
#   vertices = Vec{3,PT}((verts[1], verts[2], verts[3]))
#   tangents = Vec{2,PT}((verts[1]-verts[3], verts[2]-verts[3]))
#   normals, volume = _normals(tangents, Val{1})
#   FlatCellNM(vertices, tangents, normals, volume)
# end

# function patch{D1,U,T,D}(verts::Vec{D1,Vec{U,T}}, ::Type{Val{D}})
#     a = Array(verts)
#     patch(a, Val{2})
# end

# function patch{D1}(verts::Vec{D1})
#     patch(Array(verts), Val{D1-1})
# end



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
        r += mani.tangents[i] * u[i]
    end
    return r
end

function barytocart(mani::Array, u)
  r = mani[:,end]
  for i in 1 : size(mani,2)-1
    r += (mani[:,i] - mani[:,end]) * u[i]
  end
  return r
end

function carttobary{U,D,C,N,T}(p::FlatCellNM{U,D,C,N,T}, cart)

    G = [dot(p.tangents[i], p.tangents[j]) for i in 1:D, j in 1:D]

    w = [dot(p.tangents[i], cart - p.vertices[end]) for i in 1:D]
    u = G \ w

    return Vec(u)
end
