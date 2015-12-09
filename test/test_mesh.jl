using FixedSizeArrays
using Base.Test
MUT = CompScienceMeshes

rectangle = MUT.meshrectangle(1.0, 1.0, 1.0);

@test MUT.vertices(rectangle, Vec(1,2,3)) == [
    Point(0.0, 0.0, 0.0),
    Point(0.0, 1.0, 0.0),
    Point(1.0, 0.0, 0.0)
]

@test MUT.numvertices(rectangle) == 4
@test MUT.numcells(rectangle) == 2
@test MUT.dimension(rectangle) == 2

faces = MUT.cells(rectangle, 2)
edges = MUT.cells(rectangle, 1)
verts = MUT.cells(rectangle, 0)

@test MUT.numcells(verts) == 4
@test MUT.numcells(edges) == 5
@test MUT.numcells(faces) == 2

Λ  = MUT.connectivity(rectangle, 0, verts, edges)
Σᵀ = MUT.connectivity(rectangle, 1, edges, faces)

@test size(Λ)  == (5,4)
@test size(Σᵀ) == (2,5)

@test_approx_eq(norm(Σᵀ*Λ, Inf), 0)

bnd = MUT.boundary(rectangle)

@test MUT.dimension(bnd) == 1
@test MUT.numvertices(bnd) == 4
@test MUT.numcells(bnd) == 4

vtoe, num_vtoe = MUT.vertextocellmap(edges)
vtof, num_vtof = MUT.vertextocellmap(faces)

@test size(vtoe) == (4,3)
@test size(vtof) == (4,2)

@test size(num_vtoe) == (4,)
@test size(num_vtof) == (4,)

file = joinpath(Pkg.dir("CompScienceMeshes"),"test","sphere.in")
sphere = MUT.meshfromfile(file)
@test MUT.numvertices(sphere) == 335
@test MUT.numcells(sphere) == 666
