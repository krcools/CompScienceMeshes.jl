using CompScienceMeshes

m1 = meshrectangle(1.0, 0.5, 0.1, 3)
m1 = CompScienceMeshes.translate(m1, [0.0, 0.5, 0.0])
m2 = meshrectangle(1.0, 0.5, 0.1, 3)

nc1 = numcells(m1)
nc2 = numcells(m2)

nv1 = numvertices(m1)
nv2 = numvertices(m2)

M1 = weld(m1,m2)
@assert numcells(M1) == nc1 + nc2
@assert numvertices(M1) == nv1 + nv2 - 11

M2 = weld(m1,-m2)
@assert numcells(M2) == nc1 + nc2
@assert numvertices(M2) == nv1 + nv2 - 11

@assert CompScienceMeshes.isoriented(M1)
@assert !CompScienceMeshes.isoriented(M2)

e1 = skeleton(m1,1)
e2 = skeleton(m2,1)
E1 = skeleton(M1,1)

S11 = CompScienceMeshes.embedding(e1,E1)
S21 = CompScienceMeshes.embedding(e2,E1)
