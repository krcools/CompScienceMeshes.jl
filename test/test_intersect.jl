using Base.Test
using FixedSizeArrays

using CompScienceMeshes
using BoundaryElements

const e0 = Point(0.0,0.0,0.0)
const e1 = Point(1.0,0.0,0.0)
const e2 = Point(0.0,1.0,0.0)
const e3 = Point(0.0,0.0,0.1)

p = patch([e0, e1], Val{1})
q = patch([(e0+e1)/2, (-e0+3*e1)/2], Val{1})

pq = intersection(p,q)
@test volume(pq) == 1/2
@test pq.vertices[1] == q.vertices[1]
@test pq.vertices[2] == p.vertices[2]
@test 2*pq.tangents[1] == p.tangents[1] == q.tangents[1]
