using Test
using CompScienceMeshes

verts = [
    point(0,0,0),
    point(1,0,0),
    point(0,1,0)]

nodes = [
    (1,),
    (2,),
    (3,),]

edges = [
    (2,1,),
    (2,3,),
    (3,1,),]

faces = [
    (-1,2,3,),]

T = eltype(verts[1])
P = eltype(verts)
m = CompScienceMeshes.SComplex2D{T,P}(verts,nodes,edges,faces)

@test numcells(m) == 1

ch = chart(m, first(cells(m)))
@test ch.vertices[1] == verts[1]
@test ch.vertices[2] == verts[2]
@test ch.vertices[3] == verts[3]
