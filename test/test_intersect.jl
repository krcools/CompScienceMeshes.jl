using Base.Test
using FixedSizeArrays

using CompScienceMeshes

T = Float64
e0, e1, e2, e3 = euclidianbasis(T, 3)

p = patch([e0, e1], Val{1})
q = patch([(e0+e1)/2, (-e0+3*e1)/2], Val{1})

pq = intersection(p,q)
@test length(pq) == 1
pq = pq[1]
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

# This case caused NaNs in the RWG BC Gram matrix. Turned out to be caused
# by computing the intersection between two parallel lines. This needs a
# more robust solution in the long run.
p = Vec(-0.819550109062379,-0.394674531365384,-0.415415013001886)
q = Vec(-0.66680743687093,-0.618706884644307,-0.415415013001886)
r = Vec(-0.832690784854388,-0.535137873473586,-0.142314838273285)

a = Vec(-0.66680743687093,-0.618706884644307,-0.415415013001886)
b = Vec(-0.749749110862659,-0.5769223790589465,-0.27886492563758547)
c = Vec(-0.7730161102625658,-0.5161730964944257,-0.32438162142568566)

isct = sutherlandhodgman([p,q,r],[a,b,c])

@assert isct == [
    Vec(-0.7730161102625658,-0.5161730964944257,-0.32438162142568566)
    Vec(-0.66680743687093,-0.618706884644307,-0.415415013001886)
    Vec(-0.66680743687093,-0.618706884644307,-0.415415013001886)
    Vec(-0.749749110862659,-0.5769223790589465,-0.27886492563758547)]
