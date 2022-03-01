using Test
using CompScienceMeshes

for T in [Float32, Float64]
    mesh = meshrectangle(T(1.0), T(1.0), T(1.0))
    faces = skeleton(mesh, 2)
    local verts = vertices(mesh, faces.faces[1])
    p = simplex(verts)

    @test p.vertices == [
        point(T, 0.0, 0.0, 0.0),
        point(T, 0.0, 1.0, 0.0),
        point(T, 1.0, 0.0, 0.0)]
    @test p.tangents == [
        point(T, -1.0, 0.0, 0.0),
        point(T, -1.0, 1.0, 0.0)]
    @test p.normals == [
        point(T, 0.0, 0.0, -1.0)]
    @test p.volume == T(0.5)
end
