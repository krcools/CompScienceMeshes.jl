module CompScienceMeshes

using FixedSizeArrays

import Base.getindex

include("fsa_extensions.jl")

include("mesh.jl")
include("gmsh.jl")
include("primitives.jl")
include("simplex.jl")
include("geometry.jl")
include("patches.jl")
include("overlap.jl")
include("submesh.jl")
include("meshpoints.jl")
include("baryref.jl")
include("intersect.jl")
include("isinside.jl")
include("weld.jl")

include("utils.jl")

end # module
