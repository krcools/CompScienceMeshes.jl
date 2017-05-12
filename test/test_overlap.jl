using Base.Test
using CompScienceMeshes

p = simplex(
    point(0,0,0),
    point(1,0,0),
    point(0,1,0))

q1 = simplex(
    point(0.6, 0.6, 0.0),
    point(1.6, 0.6, 0.0),
    point(0.6, 1.6, 0.0),
)

q2 = simplex(
    point(0.4, 0.4, 0.0),
    point(1.4, 0.4, 0.0),
    point(0.4, 1.4, 0.0),
)

@test overlap(p, q1) == false
@test overlap(p, q2) == true

q1 = simplex(
    point(0.949726718617726,-0.278864925637586,0.142314838273285),
    point(0.989821441880933,0.0,-0.142314838273285),
    point(0.989821441880933,0.0,0.142314838273285),
)

q2 = simplex(
    point(0.949726718617726,-0.278864925637586,-0.142314838273285),
    point(0.989821441880933,0.0,-0.142314838273285),
    point(0.949726718617726,-0.278864925637586,0.142314838273285),
)

@test overlap(q1, q2) == false

## make sure the submesh function work for 1D meshes


l1 = meshsegment(1.0,1/2)
vt = skeleton(l1,0)
bd = boundary(l1)

overlaps = overlap_gpredicate(bd)
pred1 = c -> overlaps(simplex(vertices(vt,c)))
@test pred1(SVector(1))
@test !pred1(SVector(2))
@test pred1(SVector(3))


# test a case where the segments are:
#   not of unit length
#   colinear and opposite
#   meet in a common point
ch1 = simplex(point(1/3,0,0), point(1/3,1/3,0))
ch2 = simplex(point(1/3,1/3,0), point(1/3,2/3,0))
@test !overlap(ch1, ch2)
