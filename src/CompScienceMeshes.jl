module CompScienceMeshes

using FixedSizeArrays

import Base.getindex

#include("fsa_extensions.jl")

export euclidianbasis, defaultpointtype, Pt, point


typealias Pt{N,T} FixedSizeArrays.Vec{N,T}

"""
    point(xs...)

Create point object of `CompScienceMeshes` default point class
and with coordinate type `Float64`
"""
point(xs...) = point(Float64, xs...)

"""
    point(type, xs...)

Create point object of `CompScienceMeshes` default point class
and with coordinate type `type`
"""
@generated function point(T::Type,xs...)
    D = length(xs)
    xp = :(Vec{$D,T}())
    for d in 1:D
        push!(xp.args, :(xs[$d]))
    end
    return xp
end

export index
index(is...) = Vec{length(is),Int}(is...)


"""
  defaultpointtype(T, dim) = Vec{T,dim}

Returns the default point type used by package `CompScienceMeshes`
"""
defaultpointtype(T, dim) = Vec{dim,T}



"""
  euclidian_basis(type, dim)

Returns an arrays of length (dim+1) containing the origin and the dim
Euclidian unit vectors with coordinate type `type`.
"""
function euclidianbasis(T::DataType, dim)
  P = Vec{dim,T}
  id = eye(dim)
  r = P[ P(id[:,i]...) for i in 1:dim ]
  z = P(zeros(T,dim)...)
  return [z; r]
end

include("combinatorics.jl")

# simplices and related algorithms
include("patches.jl")
include("overlap.jl")
include("intersect.jl")
include("isinside.jl")
include("meshpoints.jl")

# mesh component
include("mesh.jl")
include("gmsh.jl")
include("primitives.jl")
include("submesh.jl")
include("baryref.jl")
include("weld.jl")

# geometry API
include("geometry.jl")



include("utils.jl")

end # module
