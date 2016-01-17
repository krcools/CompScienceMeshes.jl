module CompScienceMeshes

using FixedSizeArrays

include("fsa_extensions.jl")

include("simplex.jl")
include("mesh.jl")
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
