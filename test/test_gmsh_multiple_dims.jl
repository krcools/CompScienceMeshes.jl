using Test
using CompScienceMeshes

fn = joinpath(@__DIR__, "assets", "world.msh")
border = CompScienceMeshes.read_gmsh_mesh(fn, physical="WorldsEnd", dimension=1)
coast  = CompScienceMeshes.read_gmsh_mesh(fn, physical="Coast", dimension=1)
sea    = CompScienceMeshes.read_gmsh_mesh(fn, physical="Sea", dimension=2)

border_vertices = skeleton(border, 0)
coast_vertices  = skeleton(coast, 0)
sea_vertices    = skeleton(sea, 0)

# interior_vertices = submesh(v -> !(v in border_vertices || v in coast_vertices), sea_vertices)
interior_vertices = submesh(sea_vertices) do mesh, i
    v = CompScienceMeshes.indices(mesh, i)
    v in cells(border_vertices) && return false
    v in cells(coast_vertices) && return false
    return true
end

@test length(border_vertices) + length(coast_vertices) + length(interior_vertices) == length(sea_vertices)
