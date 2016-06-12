using Base.Test
using FixedSizeArrays

using CompScienceMeshes

T = Float64
e0, e1, e2, e3 = euclidianbasis(T, 3)

p = patch([e0, e1], Val{1})
q = patch([(e0+e1)/2, (-e0+3*e1)/2], Val{1})

pq = intersection(p,q)
@test volume(pq) == 1/2
@test pq.vertices[1] == q.vertices[1]
@test pq.vertices[2] == p.vertices[2]
@test 2*pq.tangents[1] == p.tangents[1] == q.tangents[1]

c = (e0+e1+e2)/3
p = patch([e0,e1,e2], Val{2})
q = patch([e0,e1,c], Val{2})

qq = intersection(q,q)
pq = intersection(p,q)
qp = intersection(q,p)

# For these case where points of the subject are on the clipper or
# the other way around, the precise computationof the intersection
# is subject to numerical rounding errors. For efficiency reasons
# this something that needs fixing. It should not affect the computation
# of the interaction elements, however.
@test_approx_eq(sum([volume(r) for r in qq]), volume(q))
@test_approx_eq(sum([volume(r) for r in pq]), volume(q))
@test_approx_eq(sum([volume(r) for r in qp]), volume(q))
