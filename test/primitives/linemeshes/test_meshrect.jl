using CompScienceMeshes
using StaticArrays
using Test

@testset "Test meshrect function" begin
    m = meshrect(1.0, 1.0, 1.0)

    @test numvertices(m) == 4
    @test m.vertices[3] == SVector(1.0, 1.0)
    @test m.faces[1] == CompScienceMeshes.SimplexGraph(1, 2)
end