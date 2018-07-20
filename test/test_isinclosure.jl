using CompScienceMeshes
using Test

m = meshsegment(1.0, 1.0, 3)

p1 = point(0.5, 0, 0)
p2 = point(0.5, 1e-4, 0)
p3 = point(0.0, 0, 0)
p4 = point(-1e-4, 0, 0)

f = inclosure_gpredicate(m)

##
@test f(p1) == true
@test f(p2) == false
@test f(p3) == true
@test f(p4) == false
