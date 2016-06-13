using Base.Test
using FixedSizeArrays

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
