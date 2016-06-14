using Base.Test
using CompScienceMeshes

mesh = meshrectangle(1.0, 1.0, 1.0)
faces = skeleton(mesh, 2)
verts = vertices(mesh, faces.faces[1])
p = simplex(verts)

@test p.vertices == [
    point(0.0, 0.0, 0.0),
    point(0.0, 1.0, 0.0),
    point(1.0, 0.0, 0.0)]
@test p.tangents == [
    point(-1.0, 0.0, 0.0),
    point(-1.0, 1.0, 0.0)]
@test p.normals == [
    point(0.0, 0.0, -1.0)]
@test p.volume == 0.5
