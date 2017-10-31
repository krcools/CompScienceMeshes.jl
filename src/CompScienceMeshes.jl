module CompScienceMeshes

using StaticArrays
using Compat
import Base.getindex

export getcommonedge

# defaults
export index
export euclidianbasis, point

# default mesh creation
export mesh, readmesh, writemesh
export meshsegment, meshrectangle, meshcircle, meshsphere, meshcuboid
export subdMesh
# mesh interface
export dimension, universedimension, vertextype, coordtype
export numvertices, vertices
export numcells, cells
export boundary, skeleton, vertextocellmap
export connectivity, cellpairs # marked for deprecation

# mesh transforms
export translate, translate!, rotate, rotate!
export fliporientation!, fliporientation
export weld

# mesh refinement
export barycentric_refinement, bisecting_refinement

export  Loop_subdivision

# submesh selection and cell retrieval
export GSubdMesh,submesh
export octree, boundingbox
export overlap_gpredicate, interior_tpredicate, inclosure_gpredicate

export mapper, restriction

# link to the chart concept
export chart

export domain # marked for deprecation
export dimension, universedimension
export volume
export neighborhood
export quadpoints

export isinside, isinclosure, overlap

# specific to simplicial charts
export simplex, center, vertices
export barytocart, carttobary

export intersection

# the neighborhood concept
export barycentric

#export paramtype
export cartesian, parametric
export jacobian, tangents, normal


export trgauss, sqgauss, legendre

# marked for deprecation
export SegmentedAxis
export minmaxdist, rings, ring

@compat Pt{N,T} = StaticArrays.SVector{N,T}

include("defaults.jl")
include("utils.jl")
include("combinatorics.jl")

# quadrature rules for segements, triangles, and squares
include("quadrature/SegmentGauss.jl")
include("quadrature/TriangleGauss.jl")
include("quadrature/SquareGauss.jl")

# simplices and related algorithms
include("rectangle.jl")
include("charts.jl")
include("sphere.jl")
include("overlap.jl")
include("intersect.jl")
include("isinside.jl")
include("neighborhood.jl")
include("quadpoints.jl")

# mesh component
include("mesh.jl")
include("flipped.jl")
include("timeaxis.jl")
include("subdMesh.jl")
include("subd_shape.jl")
include("gaussquarature.jl")
include("fileio/TRI_mesh.jl")
include("fileio/readmesh.jl")
include("fileio/gmsh.jl")
include("fileio/gid.jl")
include("primitives.jl")
#include("../examples/waveguide_with_post.jl")
include("submesh.jl")
include("baryref.jl")
include("subdivision.jl")
include("weld.jl")

include("mapper.jl")
include("restrict.jl")



end # module
