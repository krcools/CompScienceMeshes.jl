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

    nodes = CSM.skeleton_fast(mesh, 0)
    @test length(nodes) == 4
end


@testitem "curvy mesh: connectivity [p=2]" begin
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

    L = connectivity(nodes, edges, identity)
    S = connectivity(mesh, edges, identity)

    @test L == [
        0 1 -2 0;
        -2 0 1 0;
        1 -2 0 0;
        0 -2 0 1;
        1 0 0 -2;
    ]

    @test S == [
        1 0;
        2 0;
        3 -2;
        0 1;
        0 3;
    ]
end


@testitem "curvy mesh boundary [p=2]" begin
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
    boundary = CSM.boundary(mesh)

    @test dimension(boundary) == 1
    @test length(boundary) == 4
end


@testitem "curvy mesh bisecting refinement [p=2]" begin
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
    @show typeof(mesh)
    @test mesh isa CSM.AbstractMesh{3,3}
    fine = CSM.bisecting_refinement(mesh; order=3)
    
    @test dimension(fine) == 2
    @test length(fine) == 18
end