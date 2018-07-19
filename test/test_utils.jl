using Test
using FixedSizeArrays
CM = CompScienceMeshes

b = CM.sortneighbors([4,3,2,1], (x,y)-> y == x+1)
@test b == [1,2,3,4]

b = CM.sortneighbors([1,2,3,4], (x,y)->y == x+1)
@test b == [1,2,3,4]

b = CM.sortneighbors([4,5,6,1,2,3], (x,y)->y == x+1)
@test b == [1,2,3,4,5,6]

b = CM.sortneighbors([4,5,6,3,2,1], (x,y)->y == x+1)
@test b == [1,2,3,4,5,6]

b = CM.sortneighbors(randperm(10), (x,y)-> y == x+1)
@test b == collect(1:10)

a = [(1,2), (3,4), (2,3), (0,1)]
pred = (x,y) -> x[2] == y[1] || x[1] == y[2]
b = CM.sortneighbors(a, pred)
@test b == [(0,1), (1,2), (2,3), (3,4)]

mesh = CM.meshrectangle(1.0, 1.0, 0.5)
vtof, vton = CM.vertextocellmap(mesh)
supp = mesh.faces[vec(vtof[5,1:vton[5]])]
commonedge(x,y) = length(intersect(x,y)) == 2
supp = CM.sortneighbors(supp, commonedge)
@test length(supp) == 6
@test supp == Vec{3,Int}[
    Vec(2,5,4),
    Vec(2,3,5),
    Vec(3,6,5),
    Vec(5,6,8),
    Vec(5,8,7),
    Vec(4,5,7)
]
@test CM.isclosed(supp, commonedge) == true
