using Base.Test

using CompScienceMeshes

e0 = point(0,0,0)
e1 = point(1,0,0)
e2 = point(0,1,0)
e3 = point(0,0,1)

p = simplex(e0, e1)
q = simplex((e0+e1)/2, (-e0+3*e1)/2)

pq = intersection(p,q)
@test length(pq) == 1
pq = pq[1]
@test volume(pq) == 1/2
@test pq.vertices[1] == q.vertices[1]
@test pq.vertices[2] == p.vertices[2]
@test 2*pq.tangents[1] == p.tangents[1] == q.tangents[1]

c = (e0+e1+e2)/3
p = simplex(e0, e1, e2)
q = simplex(e0, e1, c)

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
p = point(-0.819550109062379,-0.394674531365384,-0.415415013001886)
q = point(-0.66680743687093,-0.618706884644307,-0.415415013001886)
r = point(-0.832690784854388,-0.535137873473586,-0.142314838273285)

a = point(-0.66680743687093,-0.618706884644307,-0.415415013001886)
b = point(-0.749749110862659,-0.5769223790589465,-0.27886492563758547)
c = point(-0.7730161102625658,-0.5161730964944257,-0.32438162142568566)

isct = CompScienceMeshes.sutherlandhodgman([p,q,r],[a,b,c])

@assert isct == [
    point(-0.7730161102625658,-0.5161730964944257,-0.32438162142568566)
    point(-0.66680743687093,-0.618706884644307,-0.415415013001886)
    point(-0.66680743687093,-0.618706884644307,-0.415415013001886)
    point(-0.749749110862659,-0.5769223790589465,-0.27886492563758547)]
