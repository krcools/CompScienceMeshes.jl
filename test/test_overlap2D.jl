using Test
using CompScienceMeshes
using StaticArrays

p = simplex(
    point(0,0),
    point(1,0),
    point(0,1))

q1 = simplex(
    point(0.6, 0.6),
    point(1.6, 0.6),
    point(0.6, 1.6),
)

q2 = simplex(
    point(0.4, 0.4),
    point(1.4, 0.4),
    point(0.4, 1.4),
)

@test overlap(p, q1) == false
@test overlap(p, q2) == true


# test a case where the segments are:
#   not of unit length
#   colinear and opposite
#   meet in a common point
ch1 = simplex(point(1/3,0), point(1/3,1/3))
ch2 = simplex(point(1/3,1/3), point(1/3,2/3))
@test !overlap(ch1, ch2)