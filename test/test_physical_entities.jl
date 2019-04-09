using Test
using CompScienceMeshes

fn = joinpath(@__DIR__, "assets", "thick_cladding.msh")
G01 = CompScienceMeshes.read_gmsh_mesh(fn, physical="Gamma01")
G02 = CompScienceMeshes.read_gmsh_mesh(fn, physical="Gamma02")
G12 = CompScienceMeshes.read_gmsh_mesh(fn, physical="Gamma12")

G = CompScienceMeshes.read_gmsh_mesh(fn)

@test numcells(G01) + numcells(G02) + numcells(G12) == numcells(G)
