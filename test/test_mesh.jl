CSM = CompScienceMeshes

using FixedSizeArrays
using Base.Test

rectangle = CSM.meshrectangle(1.0, 1.0, 1.0);

@test CSM.vertices(rectangle, Vec(1,2,3)) == [
    Point(0.0, 0.0, 0.0),
    Point(0.0, 1.0, 0.0),
    Point(1.0, 0.0, 0.0)
]

@test CSM.numvertices(rectangle) == 4
@test CSM.numcells(rectangle) == 2
@test CSM.dimension(rectangle) == 2

faces = CSM.cells(rectangle, 2)
edges = CSM.cells(rectangle, 1)
verts = CSM.cells(rectangle, 0)

@test length(verts) == 4
@test length(edges) == 5
@test length(faces) == 2

Λ  = CSM.connectivity(rectangle, 0, verts, edges)
Σᵀ = CSM.connectivity(rectangle, 1, edges, faces)

@test size(Λ)  == (5,4)
@test size(Σᵀ) == (2,5)

@test_approx_eq(norm(Σᵀ*Λ, Inf), 0)

bnd = CSM.boundary(rectangle)

@test CSM.dimension(bnd) == 1
@test CSM.numvertices(bnd) == 4
@test CSM.numcells(bnd) == 4

vtoe, num_vtoe = CSM.vertextocellmap(rectangle, edges)
vtof, num_vtof = CSM.vertextocellmap(rectangle, faces)

@test size(vtoe) == (4,3)
@test size(vtof) == (4,2)

@test size(num_vtoe) == (4,)
@test size(num_vtof) == (4,)
