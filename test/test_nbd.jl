using CompScienceMeshes
using Test

for T in [Float32, Float64]
    local p1 = point(T,1,0,0)
    local p2 = point(T,0,1,0)
    local p3 = point(T,0,0,1)
    local p4 = point(T,0,0,0)

    local tet = simplex(p1,p2,p3,p4)
    c = neighborhood(tet, point(T,1,1,1)/4)

    @test cartesian(c) ≈ point(T, 0.25, 0.25, 0.25)
    @test parametric(c) ≈ point(T, 0.25, 0.25, 0.25)

    @test length(cartesian(c)) == 3
    # @test c[1] ≈ 0.25
    # @test c[2] ≈ 0.25
    # @test c[3] ≈ 0.25

    @test jacobian(1.234) ≈ 1

    @test tangents(c,1) ≈ point(T,1,0,0)
    @test tangents(c,2) ≈ point(T,0,1,0)
    @test tangents(c,3) ≈ point(T,0,0,1)

    tri = simplex(p1,p2,p4)
    c_tri = neighborhood(tri, point(T,1,1)/3)

    @test normal(c_tri) ≈ point(T,0,0,1)
end