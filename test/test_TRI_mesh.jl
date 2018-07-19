using CompScienceMeshes
using Test

filename = joinpath(dirname(@__FILE__),"from_hollow_waveguide.tri.txt")
Γ = CompScienceMeshes.read_TRI_mesh(filename)

@test numvertices(Γ) == 396
@test numcells(Γ) == 300
