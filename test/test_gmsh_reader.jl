using CompScienceMeshes
using Base.Test

fn = joinpath(dirname(@__FILE__),"mesa.msh")
m = CompScienceMeshes.read_gmsh_mesh(fn)

@test numcells(m) == 1132
@test numvertices(m) == 626
@test norm(vertices(m,200) - [0.454487,-0.0982151,-0.1]) < 1e-3
