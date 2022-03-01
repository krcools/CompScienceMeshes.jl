using CompScienceMeshes
using Test

for T in [Float32, Float64]
    faces = readmesh(joinpath(@__DIR__, "assets/rectangle32.in"), T=T)
    edges = skeleton(faces,1)
    verts = skeleton(faces,0)
    bnd = boundary(faces)

    int_edges = submesh(interior_tpredicate(faces), edges)
    # pred = interior_tpredicate(int_edges)

    bnd_verts = skeleton(bnd,0)
    function pred(node)
        !(node in cells(bnd_verts))
    end
    int_verts = submesh(pred, verts)

    @test numcells(int_verts) == 9
    @test numcells(edges) - numcells(int_edges) == 16
end