using CompScienceMeshes

h = 1/6
rect = meshrectangle(1.0, 1.0, h, 3)
hole = meshrectangle(1/3, 1/3, h, 3)
translate!(hole, point(1/3, 1/3, 0))

pred = overlap_gpredicate(hole)
all_faces = submesh(c->!pred(chart(rect,c)), rect)

all_edges = skeleton(all_faces, 1)
all_verts = skeleton(all_faces, 0)

bnd_edges = boundary(all_faces)
bnd_verts = skeleton(bnd_edges, 0)

not(f) = x -> !f(x...)

onbnd1 = overlap_tpredicate(bnd_edges)
onbnd0 = overlap_tpredicate(bnd_verts)

interior_edges = submesh(not(onbnd1), all_edges)
interior_verts = submesh(not(onbnd0), all_verts)

@show numcells(interior_edges)
@show numcells(interior_verts)

D0 = connectivity(interior_verts, interior_edges)
D1 = connectivity(interior_edges, all_faces)

nullity(A) = size(A,2) - rank(A')
genus = nullity(full(D1)) - rank(full(D0))
