using Base.Test
using CompScienceMeshes

rectangle = meshrectangle(1.0, 1.0, 1.0);

@test vertices(rectangle,1) == point(0.0, 0.0, 0.0)
@test vertices(rectangle,2) == point(0.0, 1.0, 0.0)
@test vertices(rectangle,3) == point(1.0, 0.0, 0.0)

@test numvertices(rectangle) == 4
@test numcells(rectangle) == 2
@test dimension(rectangle) == 2

faces = cells(rectangle, 2)
edges = cells(rectangle, 1)
verts = cells(rectangle, 0)

@test numcells(verts) == 4
@test numcells(edges) == 5
@test numcells(faces) == 2

Λ  = connectivity(rectangle, 0, verts, edges)
Σᵀ = connectivity(rectangle, 1, edges, faces)

@test size(Λ)  == (5,4)
@test size(Σᵀ) == (2,5)

@test_approx_eq(norm(Σᵀ*Λ, Inf), 0)

bnd = boundary(rectangle)

@test dimension(bnd) == 1
@test numvertices(bnd) == 4
@test numcells(bnd) == 4

vtoe, num_vtoe = vertextocellmap(edges)
vtof, num_vtof = vertextocellmap(faces)

@test size(vtoe) == (4,3)
@test size(vtof) == (4,2)

@test size(num_vtoe) == (4,)
@test size(num_vtof) == (4,)

file = Pkg.dir("CompScienceMeshes","test","sphere.in")
sphere = meshfromfile(file)
@test numvertices(sphere) == 335
@test numcells(sphere) == 666
