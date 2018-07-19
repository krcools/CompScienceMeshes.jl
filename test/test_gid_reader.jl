using CompScienceMeshes
using Test

fn = joinpath(dirname(@__FILE__),"FPLZX.msh")
m = CompScienceMeshes.load_gid_mesh(fn)

@test numcells(m) == 896
@test numvertices(m) == 450
@test vertices(m,1) ≈ [0,0,7]
