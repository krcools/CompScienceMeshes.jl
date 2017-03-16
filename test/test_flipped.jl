using CompScienceMeshes
using Base.Test

V = [
    point(0,0,0),
    point(1,0,0),
    point(1,1,0),
    point(0,1,0)]

F = [
    index(1,2,3),
    index(1,3,4)]

m = Mesh(V,F)


f = CompScienceMeshes.FlippedMesh(m)
g = -f
@test g === m

cells(f,1)
@test cells(f,1) == index(2,1,3)
@test cells(f,2) == index(3,1,4)
