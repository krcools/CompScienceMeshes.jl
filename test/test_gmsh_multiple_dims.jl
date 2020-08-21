using Test
using CompScienceMeshes

fn = joinpath(@__DIR__, "assets", "world.msh")
border = CompScienceMeshes.read_gmsh_mesh(fn, physical="WorldsEnd", dimension=1)
coast  = CompScienceMeshes.read_gmsh_mesh(fn, physical="Coast", dimension=1)
sea    = CompScienceMeshes.read_gmsh_mesh(fn, physical="Sea", dimension=2)

border_vertices = skeleton(border, 0)P
coast_vertices  = skeleton(coast, 0)
sea_vertices    = skeleton(sea, 0)

interior_vertices = submesh(v -> !(v in border_vertices || v in coast_vertices), sea_vertices)

@test length(border_vertices) + length(coast_vertices) + length(interior_vertices) == length(sea_vertices)
