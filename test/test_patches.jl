using Base.Test
using FixedSizeArrays
MUT = CompScienceMeshes

mesh = MUT.meshrectangle(1.0, 1.0, 1.0)
faces = MUT.cells(mesh, 2)
verts = MUT.vertices(mesh, faces.faces[1])
p = MUT.patch(verts, Val{2})

@test p.vertices == Vec([
    Point(0.0, 0.0, 0.0),
    Point(0.0, 1.0, 0.0),
    Point(1.0, 0.0, 0.0)])
@test p.tangents == Vec([
    Point(-1.0, 0.0, 0.0),
    Point(-1.0, 1.0, 0.0)])
@test p.normals == Vec([
    Point(0.0, 0.0, -1.0)])
@test p.volume == 0.5
