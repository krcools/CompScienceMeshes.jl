using Base.Test
using FixedSizeArrays

import CompScienceMeshes
CM = CompScienceMeshes;

p = CM.patch(
  [
    Point(0.0, 0.0, 0.0),
    Point(1.0, 0.0, 0.0),
    Point(0.0, 1.0, 0.0),
  ], Val{2})

q1 = CM.patch(
  [
    Point(0.6, 0.6, 0.0),
    Point(1.6, 0.6, 0.0),
    Point(0.6, 1.6, 0.0),
  ], Val{2})

q2 = CM.patch(
  [
    Point(0.4, 0.4, 0.0),
    Point(1.4, 0.4, 0.0),
    Point(0.4, 1.4, 0.0),
  ], Val{2})

@test CM.overlap(p, q1) == false
@test CM.overlap(p, q2) == true

q1 = CM.patch([
    Point(0.949726718617726,-0.278864925637586,0.142314838273285),
    Point(0.989821441880933,0.0,-0.142314838273285),
    Point(0.989821441880933,0.0,0.142314838273285)], Val{2})

q2 = CM.patch([
    Point(0.949726718617726,-0.278864925637586,-0.142314838273285),
    Point(0.989821441880933,0.0,-0.142314838273285),
    Point(0.949726718617726,-0.278864925637586,0.142314838273285)], Val{2})

@test CM.overlap(q1, q2) == false
