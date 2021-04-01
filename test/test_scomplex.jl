using Test
using CompScienceMeshes
using LinearAlgebra

verts = [
    point(0,0,0),
    point(1,0,0),
    point(0,1,0),
    point(1,1,0),]

nodes = [
    (1,),
    (2,),
    (3,),
    (4,)]

edges = [
    (1,2,),
    (2,3,),
    (3,1,),
    (2,4),
    (4,3),]

faces = [
    (1,2,3,),
    (-2,4,5)]

T = eltype(verts[1])
P = eltype(verts)
m = CompScienceMeshes.SComplex2D{T,P}(verts,nodes,edges,faces)

@test numcells(m) == 2

ch = chart(m, first(cells(m)))
@test ch.vertices[1] == verts[3]
@test ch.vertices[2] == verts[1]
@test ch.vertices[3] == verts[2]
@test dot(normal(center(ch)),point(0,0,1)) > 0

x1 = iterate(cells(m))
@test x1 != nothing
x2 = iterate(cells(m), x1[2])
@test x2 != nothing
c2 = x2[1]

ch2 = chart(m, c2)
@test ch2.vertices[1] == verts[4]
@test ch2.vertices[2] == verts[3]
@test ch2.vertices[3] == verts[2]
@test dot(normal(center(ch2)),point(0,0,1)) > 0

m0 = skeleton(m,0)
@test dimension(m0) == 0

m1 = skeleton(m,1)
@test dimension(m1) == 1

m2 = skeleton(m,2)
@test dimension(m2) == 2

D = connectivity(m1,m2)
@test sum(D,dims=1) == [1 0 1 1 1]

# b = CompScienceMeshes.SComplex1D{T,P}(verts, nodes, edges[[1,3,4,5]])
b = boundary(m2)
@test numcells(b) == 4
b0 = skeleton(b,0)
@test numcells(b0) == 4

@test inclosure_gpredicate(b).(verts) == [true, true, true, true]

mm2 = weld(m,m,seam=b)
mm1 = skeleton(mm2,1)
mm0 = skeleton(mm2,0)
@test numcells(mm2) == 2*numcells(m)
@test numcells(mm1) == 6
@test numcells(mm0) == 4

r_old = meshrectangle(1.0, 1.0, 0.5, 3)
r_new = CompScienceMeshes.SComplex2D(r_old)

@test numcells(r_new) == 8
@test numcells(skeleton(r_new,1)) == 16
@test numcells(skeleton(r_new,0)) == 9

b_new = boundary(r_new)
@test dimension(b_new) == 1
@test numcells(b_new) == 8

# r_new1 = skeleton(r_new,1)
# p = interior_tpredicate(r_new)
# for e in cells(r_new1)
#     println(cartesian(center(chart(r_new1,e))))
#     println(p(e))
# end

rr = weld(r_new, r_new, seam=boundary(r_new))
rr1 = skeleton(rr,1)
rr0 = skeleton(rr,0)
@test numcells(rr) == 16
@test numcells(rr1) == 24
@test numcells(rr0) == 10

fm = -m
C = collect(cells(fm))
ch1 = chart(fm,C[1])
ch2 = chart(fm,C[2])

@test dot(normal(center(ch1)), point(0,0,1)) < 0
@test dot(normal(center(ch2)), point(0,0,1)) < 0

dl = weld(m,-m,seam=boundary(m))
@test numcells(dl) == 4
@test numcells(skeleton(dl,1)) == 6
@test numcells(skeleton(dl,0)) == 4

#Test interior
i_old = CompScienceMeshes.interior(r_old)
i_new = CompScienceMeshes.interior(r_new)
@test dimension(i_new) == 1
@test numcells(i_new) == 8
@test numcells(i_new) == numcells(i_old)
@test numcells(i_new) + numcells(b_new) == numcells(skeleton(r_new,1))

sn = skeleton(i_new,0)
so = skeleton(i_old,0)
@assert numcells(sn) == numcells(so)
for (i,v) in enumerate(cells(sn))
    @test v[1] == cells(so)[i][1]
end

#Test Mesh(SComplex)
c2 = meshsphere(1.0,0.8)
c1 = skeleton(c2,1)
c0 = skeleton(c2,0)

d2 = CompScienceMeshes.SComplex2D(c2)
d1 = skeleton(d2,1)
d0 = skeleton(d2,0)

e2 = Mesh(d2)
e1 = Mesh(d1)
e0 = Mesh(d0)

@test e2.vertices == c2.vertices
@test e1.vertices == c1.vertices
@test e0.vertices == c0.vertices

@test e2.faces == c2.faces
@test e1.faces == c1.faces
@test e0.faces == c0.faces

f2 = CompScienceMeshes.SComplex2D(e2)

@test f2.faces == d2.faces
@test f2.edges == d2.edges
@test f2.nodes == d2.nodes
