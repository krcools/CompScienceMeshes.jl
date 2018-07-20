using CompScienceMeshes
using Test

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

T = [index(2,1,3), index(3,1,4)]
for (i,_c) in enumerate(cells(f))
    @test _c == T[i]
end

# @test f.faces[1] == index(2,1,3)
# @test f.faces[2] == index(3,1,4)
