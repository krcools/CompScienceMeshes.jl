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

## zero-dimensional cells: the closure of a point is the point itself
b = boundary(meshsegment(1.0,1/4,3))
g = inclosure_gpredicate(b)

@test g(point(0.0,0,0)) == true #start point
@test g(point(1.0,0,0)) == true # end point
@test g(point(0.5,0,0)) == false # interior point
@test g(point(0.0,1.0,0)) == false # point off the segment
