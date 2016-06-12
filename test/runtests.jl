module PkgTests
import CompScienceMeshes

using FixedSizeArrays
using Base.Test

@test CompScienceMeshes.point(1,2,3) === Vec{3,Float64}(1.0,2.0,3.0)
@test CompScienceMeshes.point(Int,1,2,3) == Vec{3,Int}(1,2,3)

include("test_mesh.jl")
include("test_geometry.jl")
include("test_patches.jl")
include("test_submesh.jl")
include("test_overlap.jl")
include("test_intersect.jl")
include("test_sh_intersection.jl")
include("test_isinside.jl")

end
