module CompScienceMeshes

using StaticArrays
using Compat
import Base.getindex

# defaults
export index
export euclidianbasis, point

# default mesh creation
export mesh, readmesh, writemesh
export meshsegment, meshrectangle, meshcircle, meshsphere, meshcuboid

# mesh interface
export dimension, universedimension, vertextype, coordtype
export numvertices, vertices
export numcells, cells
export boundary, skeleton, vertextocellmap
export connectivity, cellpairs # marked for deprecation

# mesh transforms
export translate, translate!, rotate, rotate!
export fliporientation!, fliporientation

# mesh refinement
export barycentric_refinement, bisecting_refinement

# link to the chart concept
export chart


export domain # marked for deprecation
export dimension, universedimension
export tangents, volume
export cartesian, jacobian, neighborhood
export intersection
export isinside, isinclosure, overlap

# specific to simplicial charts
export simplex, center, vertices
export barycentric, barytocart, carttobary

# submesh selection and cell retrieval
export submesh, octree, boundingbox
export overlap_gpredicate, interior_tpredicate, inclosure_gpredicate

export neighborhood
export paramtype
export cartesian, parametric
export jacobian, tangents, normal

export quadpoints
export trgauss, sqgauss, legendre

# marked for deprecation
export SegmentedAxis
export minmaxdist, rings, ring

#typealias Pt{N,T} StaticArrays.SVector{N,T}
@compat Pt{N,T} = StaticArrays.SVector{N,T}

include("defaults.jl")
include("utils.jl")

# quadrature rules for segements, triangles, and squares
include("quadrature/SegmentGauss.jl")
include("quadrature/TriangleGauss.jl")
include("quadrature/SquareGauss.jl")

# simplices and related algorithms
include("charts.jl")
include("overlap.jl")
include("intersect.jl")
include("isinside.jl")
include("neighborhood.jl")
include("quadpoints.jl")

# mesh component
include("mesh.jl")
include("flipped.jl")
include("timeaxis.jl")

include("gmsh.jl")
include("gid.jl")
include("primitives.jl")
include("submesh.jl")
include("baryref.jl")
include("weld.jl")

end # module
