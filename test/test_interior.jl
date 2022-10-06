using CompScienceMeshes
using Test

T = Float64
for T in [Float32, Float64]
    local fn = joinpath(@__DIR__, "assets/rectangle32.in")
    local faces = readmesh(fn, T=T)
    local edges = skeleton(faces,1)
    local verts = skeleton(faces,0)
    local bnd = boundary(faces)

    local int_edges = submesh(interior_tpredicate(faces), edges)
    local bnd_verts = skeleton(bnd,0)

    local cells_bnd_verts = [CompScienceMeshes.indices(bnd_verts,i) for i in bnd_verts]
    function pred(mesh, node)
        inds = CompScienceMeshes.indices(mesh, node)
        !(inds in cells_bnd_verts)
    end
    local int_verts = submesh(pred, verts)

    @test numcells(int_verts) == 9
    @test numcells(edges) - numcells(int_edges) == 16
end