using Test
using LinearAlgebra
using CompScienceMeshes

@testset "meshbunny" begin
    h = 0.1
    m = meshbunny(h)
    @test m isa CompScienceMeshes.Mesh
    @test CompScienceMeshes.isoriented(m)
    nv, ne, nf = length(skeleton(m, 0)), length(skeleton(m, 1)), numcells(m)
    @test nv - ne + nf == 2                     # closed and of genus 0
    @test length(boundary(m)) == 0
    @test maximum(reduce((a, b) -> max.(a, b), vertices(m)) - reduce((a, b) -> min.(a, b), vertices(m))) ≈ 1 rtol=0.05
    # outward orientation: positive enclosed volume
    @test sum(dot(c[1], cross(c[2], c[3])) for c in (vertices(m)[i.indices] for i in cells(m))) > 0

    m3 = meshbunny(h; order=3)
    @test m3 isa CompScienceMeshes.CurvilinearMesh
    @test meshorder(m3) == 3
    @test numcells(m3) == numcells(m)
    @test CompScienceMeshes.isoriented(m3)
    # the cubic patches are not folded: their normals stay within 90° of the flat normals
    pts = [point(u, v) for u in 0:0.25:1 for v in 0:0.25:1 if u + v <= 1]
    @test all(1:numcells(m3)) do i
        nf = normal(center(chart(m, cells(m)[i])))
        all(dot(normal(neighborhood(chart(m3, i), p)), nf) > 0 for p in pts)
    end
end
