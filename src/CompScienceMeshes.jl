module CompScienceMeshes

using StaticArrays
import Base.getindex

export getcommonedge

# defaults
export index
export euclidianbasis, point

export mesh, readmesh, writemesh
export dimension, universedimension, vertextype, coordtype
export numvertices, vertices
export numcells, cells
export translate, translate!, rotate, rotate!
export fliporientation!, fliporientation
export boundary, skeleton, vertextocellmap, connectivity, cellpairs
export barycentric_refinement, bisecting_refinement
export chart

# primitives
export meshsegment, meshrectangle, meshcircle, meshsphere, meshcuboid
#export meshwaveguidepost  # Deprecate

export domain
export simplex
export dimension, universedimension
export vertextype
export vertices, tangents, volume
export barytocart, carttobary, centroid
export cartesian, jacobian, neighborhood
export intersection
export isinside, isinclosure, overlap

export submesh, octree, boundingbox
export overlap_gpredicate, interior_tpredicate, inclosure_gpredicate

export neighborhood
export meshpoints
export paramtype
export cartesian, parametric, barycentric
export jacobian, tangents, utangents, normal

export quadpoints
export trgauss, sqgauss, legendre

# TODO: remove this and use default mesh interface
export SegmentedAxis
export minmaxdist, rings, ring

typealias Pt{N,T} StaticArrays.SVector{N,T}

include("defaults.jl")
include("combinatorics.jl")

# quadrature rules for segements, triangles, and squares
include("quadrature/SegmentGauss.jl")
include("quadrature/TriangleGauss.jl")
include("quadrature/SquareGauss.jl")

# simplices and related algorithms
include("patches.jl")
include("overlap.jl")
include("intersect.jl")
include("isinside.jl")
include("meshpoints.jl")
include("quadpoints.jl")

# mesh component
include("mesh.jl")
include("flipped.jl")
include("timeaxis.jl")

include("gmsh.jl")
include("gid.jl")
include("primitives.jl")
include("../examples/waveguide_with_post.jl")
include("submesh.jl")
include("baryref.jl")
include("weld.jl")

# geometry API
include("geometry.jl")



include("utils.jl")

end # module
