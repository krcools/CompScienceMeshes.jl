# overload, don't overwrite!
import CompScienceMeshes.dimension

# Factory methods
export patch

# Interface for any Patch
export barytocart, vertices, dimension, volume

# Interface for FlatPatch
export tangents, dimension, centroid


abstract Patch{T}
abstract FlatCell{T} <: Patch{T}

vertices(p::FlatCell) = p.vertices
tangents(p::FlatCell) = p.tangents
volume(p::FlatCell) = p.volume
dimension(p::FlatCell) = size(p.vertices, 2) - 1


# U: the dimension of the universe
# D: the dimension of the manifold
# N: the number of vertices
# T: the type of the coordinates
# C: the complimentary dimension (should always be U-D)
immutable FlatCellNM{U,D,C,N,T} <: FlatCell{T}
    vertices::Vec{N,Point{U,T}}
    tangents::Vec{D,Point{U,T}}
    normals::Vec{C,Point{U,T}}
    volume::T
end

dimension{U,D,C,N,T}(p::FlatCellNM{U,D,C,N,T}) = D

function patch{U,D,T}(verts::Array{Point{U,T}, 1}, ::Type{Val{D}})

    # compute the tangents
    PT = eltype(verts)
    vertices = Vec{D+1,PT}(verts)
    tangents = Vec{D,PT}([verts[i]-verts[end] for i in 1:D])

    normals, volume = _normals(tangents, Val{U-D})

    return FlatCellNM{U,D,U-D,D+1,T}(vertices, tangents, normals, volume)
end

function patch{T}(verts::Array{Point{3,T},1}, ::Type{Val{2}})
    PT = eltype(verts)
    vertices = Vec{3,PT}(verts)
    tangents = Vec{2,PT}((verts[1]-verts[3], verts[2]-verts[3]))
    normals, volume = _normals(tangents, Val{1})
    return FlatCellNM{3,2,1,3,T}(vertices, tangents, normals, volume)
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
    T  = eltype(PT)

    metric = T[dot(tangents[i], tangents[j]) for i in 1:D, j in 1:D]
    temp = [tangents[j][i] for i in 1:U, j in 1:D]
    normals = nullspace(temp * inv(metric) * temp')
    for i in 1:C
        normals[:,i] /= norm(normals[:,i])
    end
    normals = Vec{C,PT}([Point(normals[:,i]) for i in 1:C])
    volume = sqrt(abs(det(metric))) /  D

    return normals, volume
end



function patch{U,D1,T}(verts::Vec{D1,Point{U,T}})
    patch(Array(verts), Val{D1-1})
end



function barytocart(mani::FlatCellNM, u)
    r = last(mani.vertices)
    for i in 1 : length(u)
        r += mani.tangents[i] * u[i]
    end
    return r
end
