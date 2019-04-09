using Test
using CompScienceMeshes

# G0 is the interior of the metal (no simulation required here)
# G1 is the unbounded exterior
# G2 is the bounded interior

fn = joinpath(@__DIR__, "assets", "thick_cladding.msh")
G01 = CompScienceMeshes.read_gmsh_mesh(fn, physical="Gamma01")
G02 = CompScienceMeshes.read_gmsh_mesh(fn, physical="Gamma02")
G12 = CompScienceMeshes.read_gmsh_mesh(fn, physical="Gamma12")

G = CompScienceMeshes.read_gmsh_mesh(fn)

@test numcells(G01) + numcells(G02) + numcells(G12) == numcells(G)

# Build region boundaries with normals pointing inward:
G1 = weld(-G01, -G12)
@test CompScienceMeshes.isoriented(G1)

G2 = weld(-G02, G12)
@test CompScienceMeshes.isoriented(G2)

isclosed(m) = (numcells(boundary(m)) == 0)

@test isclosed(G1)
@test isclosed(G2)
