module PkgTests

using LinearAlgebra
using InteractiveUtils

using SparseArrays
using StaticArrays
using Test

import CompScienceMeshes

@test CompScienceMeshes.point(1,2,3) === SVector{3,Float64}(1.0,2.0,3.0)
@test CompScienceMeshes.point(Int,1,2,3) == SVector{3,Int}(1,2,3)

include("test_mesh.jl")
include("test_scomplex.jl")
include("test_flipped.jl")
include("test_flipped_welding.jl")
include("test_weld.jl")
include("test_baryref.jl")
include("test_gid_reader.jl")
include("test_gmsh_reader.jl")
include("test_physical_entities.jl")

include("test_geometry.jl")

include("test_patches.jl")
include("test_submesh.jl")
include("test_overlap.jl")
include("test_intersect.jl")
include("test_sh_intersection.jl")
include("test_isinside.jl")
include("test_jctweld.jl")
include("test_isinclosure.jl")

include("test_trgauss.jl")
include("test_chartquad.jl")
include("test_spherequad.jl")
include("test_nbd.jl")

include("test_mapper.jl")

end
