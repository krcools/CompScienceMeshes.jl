using Test

using CompScienceMeshes

e0 = point(0,0,0)
e1 = point(1,0,0)
e2 = point(0,1,0)
e3 = point(0,0,1)

p = simplex(e0,e1,e2)

@test isinside(p, (e0+e1+e2)/3) == true
@test isinside(p, point(-0.1, 0.3, 0.0)) == false
@test isinside(p, point(0.3, -0.3, 0.0)) == false
@test isinside(p, point(0.6, 0.6, 0.0)) == false
@test isinside(p, point(0.4, 0.6, 0.0)) == false


@test isinclosure(p, (e0+e1+e2)/3) == true
@test isinclosure(p, point(-0.1, 0.3, 0.0)) == false
@test isinclosure(p, point(0.3, -0.3, 0.0)) == false
@test isinclosure(p, point(0.6, 0.6, 0.0)) == false
@test isinclosure(p, point(0.4, 0.6, 0.0)) == true
