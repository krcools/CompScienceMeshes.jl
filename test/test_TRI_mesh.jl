using CompScienceMeshes
using Base.Test

filename = "from_hollow_waveguide.tri.txt"
Γ = CompScienceMeshes.read_TRI_mesh(filename)

@test numvertices(Γ) == 396
@test numcells(Γ) == 300
