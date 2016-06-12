
import Base.getindex

# Factory methods
export patch

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

coordtype(p::FlatCellNM) = eltype(eltype(p.vertices))

vertices(p::FlatCellNM) = p.vertices
tangents(p::FlatCellNM) = p.tangents
normal(p::FlatCellNM) = p.normals[1]
volume(p::FlatCellNM) = p.volume
dimension(p::FlatCellNM) = length(p.vertices) - 1
getindex(p::FlatCellNM, I) = p.vertices[I]


function patch{P<:FixedArray,D}(verts::Vector{P}, ::Type{Val{D}})
  @assert length(verts) == D+1
  U = length(P)
  T = eltype(P)
  # TODO: to this in a generated function
  vertices = Vec{D+1,P}(verts)
  tangents = Vec{D,P}([verts[i]-verts[end] for i in 1:D])
  normals, volume = _normals(tangents, Val{U-D})
  FlatCellNM(vertices, tangents, normals, volume)
end

function patch{T}(verts::Vector{Vec{3,T}}, ::Type{Val{2}})
  @assert length(verts) == 3
  PT = eltype(verts)
  vertices = Vec{3,PT}((verts[1], verts[2], verts[3]))
  tangents = Vec{2,PT}((verts[1]-verts[3], verts[2]-verts[3]))
  normals, volume = _normals(tangents, Val{1})
  FlatCellNM(vertices, tangents, normals, volume)
end

function patch{D1,U,T,D}(verts::Vec{D1,Vec{U,T}}, ::Type{Val{D}})
    a = Array(verts)
    patch(a, Val{2})
end

function patch{D1}(verts::Vec{D1})
    patch(Array(verts), Val{D1-1})
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
  barytocart(mani, u) -> p
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
