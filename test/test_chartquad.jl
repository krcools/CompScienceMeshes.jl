using CompScienceMeshes
using StaticArrays
using Base.Test

p1 = point(1.0, 0.0, 0.0)
p2 = point(10.0, 0.0, 1.0)
p3 = point(0.0, 2.0, 5.0)

ch = simplex(p1,p2,p3)

dom = domain(ch)
q, u = quadpoints(dom, 7)
@which neighborhood(ch, q[1])

p, w = CompScienceMeshes.quadpoints(ch, 7);

t1 = sum(w)
t2 = 0.5*norm(cross(p1-p3, p2-p3))

@test 2 * (t1-t2) / (t1+t2) + 1  ≈ 1 
