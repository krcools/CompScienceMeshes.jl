@testitem "chart and graph consistency [p=2]" begin
    using SparseArrays

    const CSM = CompScienceMeshes

    verts = [
        point(1.0, 0.0, 0.0),
        point(0.0, 1.0, 0.0),
        point(0.0, 0.0, 0.0),
        point(0.0, 0.5000000000013304, 0.0),
        point(0.4999999999986718, 0.0, 0.0),
        point(0.5000000000013293, 0.4999999999986707, 0.0),

        point(1.0, 1.0, 0.0),
        point(1.0, 0.5, 0.0),
        point(0.5, 1.0, 0.0),
        
    ]

    I = CSM.HOSimplexGraph{3,2}
    faces = [
        I(1, 2, 3, 6, 4, 5),
        I(1, 7, 2, 8, 9, 6),
    ]

    mesh = CSM.CurvilinearMesh(verts, faces, 2, 2)
    edges = CSM.skeleton_fast(mesh, 1)
    nodes = CSM.skeleton_fast(mesh, 0)

    L = connectivity(nodes, edges, abs)
    S = connectivity(mesh, edges, abs)

    ch1 = CSM.chart(mesh, 1)
    cl1 = CSM.cell(mesh, 1)

    nnzs = nonzeros(S)
    rows = rowvals(S)
    for k = nzrange(S, 1)[1]
        local_idx = nnzs[k]
        global_idx = rows[k]

        fc1 = CSM.faces(ch1, Val{1})[local_idx]
        fc2 = CSM.chart(edges, global_idx)
        
        x1 = CSM.center(fc1) |> cartesian
        x2 = CSM.center(fc2) |> cartesian

        @test x1 ≈ x2
    end
end