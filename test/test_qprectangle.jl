using Test
using CompScienceMeshes
using StaticArrays

o = point(0,0,0)
x = point(1,0,0)
y = point(0,1,0)

rect = CompScienceMeshes.RectangleChart(o, SVector((x,y)))
ctr = neighborhood(rect, (0.5,0.5))
@test cartesian(ctr) ≈ point(0.5, 0.5, 0.0)
@test normal(ctr) ≈ point(0,0,1)

QP = quadpoints(rect, (5,6))

@test sum(qp[2] for qp in QP) ≈ 1.0

f(r) = r[1]*r[2]
@test sum(w*f(p) for (p,w) in QP) ≈ 1/4

S(r) = point(1,1,1)
sum(w*dot(S(p),normal(p)) for (p,w) in QP)
