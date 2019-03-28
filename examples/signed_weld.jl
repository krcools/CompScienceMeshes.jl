using CompScienceMeshes


m1 = meshrectangle(1.0, 0.5, 0.1, 3)
m1 = CompScienceMeshes.translate(m, [0.0, 0.5, 0.0])
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

using MATLAB
mat"""
    figure()
    hold on
    m = Mesh($(vertexarray(M2)), $(cellarray(M2)))
    patch(m)
    [c,n] = faceNormals(m)
    quiver3x(c,n)
    true"""
