using CompScienceMeshes
using Base.Test

vrts = [
    point(0.0, 0.0, 0.0),
    point(1.0, 0.0, 0.0),
    point(0.0, 1.0, 0.0),
]

idcs = [
    index(1,2,3)
]

axis = 0.5π * [1,0,0]
G1 = Mesh(vrts, idcs)
G2 = rotate(G1, axis)
G3 = rotate(G1, 2*axis)

G = weld(G1, G2);
G = weld(G1, G2, G3)

@test numcells(G) == numcells(G1) + numcells(G2) + numcells(G3)
@test numvertices(G) == 5
