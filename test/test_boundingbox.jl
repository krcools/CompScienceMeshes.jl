using CompScienceMeshes
using Test

m = meshrectangle(1.0,1.0,1.0,3)
p = chart(m, cells(m,i))

# bouding box of an array of points
c, s = boundingbox(p.vertices)
@test isa(c, eltype(p.vertices))
@test isa(s, Number)

# bouding box of a point
c, s = boundingbox(p.vertices[1])
@test c == p.vertices[1]
@test s == 0

@testset "simplex bounding-box invariants" begin
    for T in (Float32, Float64)
        q = simplex(
            point(T, -2, 1, 3),
            point(T, 4, -1, 0),
            point(T, 1, 5, 2),
        )
        center, halfsize = boundingbox(q)
        lower = center .- halfsize
        upper = center .+ halfsize
        @test all(all(lower .<= vertex .<= upper) for vertex in q.vertices)
        @test halfsize >= zero(T)
        @test center == (minimum(q.vertices) + maximum(q.vertices)) / 2
    end

    sphere = readmesh(joinpath(@__DIR__, "assets", "sphere8.in"))
    for cell in cells(sphere)
        patch = chart(sphere, cell)
        center, halfsize = boundingbox(patch)
        @test all(all(center .- halfsize .<= vertex .<= center .+ halfsize)
            for vertex in patch.vertices)
        @test halfsize >= zero(coordtype(sphere))
    end
end
