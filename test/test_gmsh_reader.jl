using CompScienceMeshes
using Base.Test

fn = joinpath(dirname(@__FILE__),"mesa.msh")
m = CompScienceMeshes.read_gmsh_mesh(fn)

# @test numcells(m) == 896
# @test numvertices(m) == 450
# @test vertices(m,1) ≈ [0,0,7]
