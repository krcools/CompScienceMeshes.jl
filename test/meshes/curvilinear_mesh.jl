@testitem "curvy mesh: skeleton [p=2]" begin
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

    @test length(edges) == 5
end