module CompScienceMeshes

using FixedSizeArrays

include("utils.jl")

include("simplex.jl")
include("mesh.jl")
include("geometry.jl")
include("patches.jl")
include("overlap.jl")
include("submesh.jl")
include("meshpoints.jl")

end # module
