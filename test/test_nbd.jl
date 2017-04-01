using CompScienceMeshes
using Base.Test

p1 = point(1,0,0)
p2 = point(0,1,0)
p3 = point(0,0,1)
p4 = point(0,0,0)

tet = simplex(p1,p2,p3,p4)
c = neighborhood(tet, point(1,1,1)/4)

@test cartesian(c) ≈ point(0.25, 0.25, 0.25)
@test parametric(c) ≈ point(0.25, 0.25, 0.25)

@test length(c) == 3
@test c[1] ≈ 0.25
@test c[2] ≈ 0.25
@test c[3] ≈ 0.25

@test jacobian(1.234) ≈ 1

@test tangents(c,1) ≈ point(1,0,0)
@test tangents(c,2) ≈ point(0,1,0)
@test tangents(c,3) ≈ point(0,0,1)

@test utangents(c,1) ≈ point(1,0,0)
@test utangents(c,2) ≈ point(0,1,0)
@test utangents(c,3) ≈ point(0,0,1)

tri = simplex(p1,p2,p4)
c_tri = neighborhood(tri, point(1,1)/3)

@test normal(c_tri) ≈ point(0,0,1)
