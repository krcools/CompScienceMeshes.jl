using Base.Test

MUT = CompScienceMeshes


T = Float64
tol = eps(T) * 10^3
const third = one(T)/3

# test creation of a rectangle
m = MUT.meshrectangle(1.0,1.0,0.5);
@test MUT.numvertices(m) == 9
@test MUT.numcells(m) == 8

# test edge detection
edges = MUT.cells(m, 1)
@test MUT.numcells(edges) == 16

# test volume pair construction
# from a list of edges
vp = MUT.buildvolmpairs(m,edges)
internal_edges = count(x->x>0, vp[2,:])
@test internal_edges == 8

# test the relative orientation routine
@test MUT.relorientation([1,2,3],[1,2,3,4]) == -4
@test MUT.relorientation([2,3,4],[1,2,3,4]) == +1
@test MUT.relorientation([1,3,2],[1,2,3,4]) == +4
@test MUT.relorientation([2,4,3],[1,2,3,4]) == -1

v = [
    Point(0.0, 0.0, 0.0),
    Point(1.0, 0.0, 0.0),
    Point(0.0, 1.0, 0.0)]
cell = MUT.patch(v, Val{2})

A  = MUT.volume(cell)
@test abs(A - 0.5) < tol

mp = MUT.meshpoint(cell,[third, third]);
r = MUT.cartesian(mp)
@test norm(r - Point(third, third, zero(T))) < tol


# repeat the test but now on a random cell
v = [
    Point(2*rand(T,3)-1),
    Point(2*rand(T,3)-1),
    Point(2*rand(T,3)-1)]
cell = MUT.patch(v, Val{2})

mp = MUT.meshpoint(cell,[third, third]);
r = MUT.cartesian(mp)
v1 = v[1]
v2 = v[2]
v3 = v[3]

@test norm(third*(v1 + v2 + v3) - r) < tol
@test abs(dot(cross(r-v1, r-v2), r-v3)) < tol

# test that the normal is directed outwards
m = MUT.meshcircle(1.0, 2π/51)
for i in 1:MUT.numcells(m)
    p = MUT.patch(m.vertices[m.faces[i]], Val{1})
    c = (p.vertices[1] + p.vertices[2]) / 2
    @test dot(c, p.normals[1]) > 0
end
