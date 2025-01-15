using CompScienceMeshes
using Test

m1 = meshrectangle(1.0, 0.5, 0.1)
m1 = CompScienceMeshes.translate(m1, [0.0, 0.5, 0.0])
m2 = meshrectangle(1.0, 0.5, 0.1)

nc1 = numcells(m1)
nc2 = numcells(m2)

nv1 = numvertices(m1)
nv2 = numvertices(m2)

M1 = weld(m1,m2)
@test numcells(M1) == nc1 + nc2
@test numvertices(M1) == nv1 + nv2 - 11

M2 = weld(m1,-m2)
@test numcells(M2) == nc1 + nc2
@test numvertices(M2) == nv1 + nv2 - 11
