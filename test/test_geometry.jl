using Base.Test
using CompScienceMeshes


T = Float64
tol = eps(T) * 10^3
const third = one(T)/3

# test creation of a rectangle
m = meshrectangle(1.0,1.0,0.5);
@test numvertices(m) == 9
@test numcells(m) == 8

# test edge detection
edges = skeleton(m, 1)
@test numcells(edges) == 16

# test volume pair construction from a list of edges
vp = cellpairs(m,edges)
internal_edges = count(x->x>0, vp[2,:])
@test internal_edges == 8

# test the relative orientation routine
@test CompScienceMeshes.relorientation(index(1,2,3),index(1,2,3,4)) == -4
@test CompScienceMeshes.relorientation(index(2,3,4),index(1,2,3,4)) == +1
@test CompScienceMeshes.relorientation(index(1,3,2),index(1,2,3,4)) == +4
@test CompScienceMeshes.relorientation(index(2,4,3),index(1,2,3,4)) == -1

o, e1, e2, e3 = euclidianbasis(Float64, 3)
cell = simplex(o, e1, e2)

A  = volume(cell)
@test abs(A - 0.5) < tol

mp = meshpoint(cell,[third, third]);
r = cartesian(mp)
@test norm(r - Vec(third, third, zero(T))) < tol

# repeat the test but now on a random cell
#P = Pt{3,T}
#randpoint() = P(2rand(T,3)-1)
randpoint() = point(2*rand(T)-1, 2*rand(T)-1, 2*rand(T)-1)
v1 = randpoint()
v2 = randpoint()
v3 = randpoint()
cell = simplex(v1, v2, v3)

mp = meshpoint(cell,[third, third]);
r = cartesian(mp)

@test norm(third*(v1 + v2 + v3) - r) < tol
@test abs(dot(cross(r-v1, r-v2), r-v3)) < tol

# test that the normal is directed outwards
m = meshcircle(1.0, 2π/51)
for i in 1:numcells(m)
    p = simplex(cellvertices(m,i))
    c = (p.vertices[1] + p.vertices[2]) / 2
    @test dot(c, p.normals[1]) > 0
end


## test cellpairs on non-oriented meshes
using CompScienceMeshes
using Base.Test

p1 = point(1,0,0)
p2 = point(0,1,0)
p3 = point(-1,0,0)
p4 = point(0,-1,0)

i1 = index(1,2,3)
i2 = index(1,4,2)

m = Mesh([p1,p2,p3,p4], [i1,i2])
e = skeleton(m,1)

@test numcells(m) == 2
@test numcells(e) == 5

cps = cellpairs(m,e)
