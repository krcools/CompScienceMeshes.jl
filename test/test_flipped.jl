using CompScienceMeshes
using Test

for U in [Float32, Float64]
    global V = [
        point(U,0,0,0),
        point(U,1,0,0),
        point(U,1,1,0),
        point(U,0,1,0)]

    C = CompScienceMeshes.SimplexGraph{3}
    global F = [
        C(index(1,2,3)),
        C(index(1,3,4))]

    global m = Mesh(V,F)


    global f = CompScienceMeshes.FlippedMesh(m)
    global g = -f
    @test g === m

    global T = [C(index(2,1,3)), C(index(3,1,4))]
    # for (i,_c) in enumerate(cells(f))
    for (i,cell) in pairs(cells(f))
        # inds = CompScienceMeshes.indices(f,i)
        @test cell == T[i]
    end
end
# @test f.faces[1] == index(2,1,3)
# @test f.faces[2] == index(3,1,4)
