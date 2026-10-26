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


@testitem "curvy mesh bisecting refinement [rectangle, p=2]" begin
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
    @show typeof(fine)

    @test length(skeleton(mesh,0)) - length(skeleton(mesh,1)) + length(mesh) == 1
    @test length(skeleton(fine,0)) - length(skeleton(fine,1)) + length(fine) == 1
    
    @test dimension(fine) == 2
    @test length(fine) == 18
end


@testitem "curvy mesh bisecting refinement [sphere, p=2]" begin
    const CSM = CompScienceMeshes

    mesh = gmshsphere(1.0, 0.5; order=2)
    @show typeof(mesh)
    @test mesh isa CSM.AbstractMesh{3,3}
    fine = CSM.bisecting_refinement(mesh; order=3)

    @test length(skeleton(mesh,0)) - length(skeleton(mesh,1)) + length(mesh) == 2
    @test length(skeleton(fine,0)) - length(skeleton(fine,1)) + length(fine) == 2
end


@testitem "curvy mesh bisecting refinement injections [rectangle, p=2]" begin
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
    fine = CSM.bisecting_refinement(mesh; order=3)
    parent = CSM.parent(fine)

    @test length(skeleton(mesh,0)) - length(skeleton(mesh,1)) + length(mesh) == 1
    @test length(skeleton(fine,0)) - length(skeleton(fine,1)) + length(fine) == 1
    
    @show length(skeleton(fine,1)) == 33

    # @test dimension(fine) == 2
    # @test length(fine) == 18

    # c1 = cells(parent)[1]
    for t in parent
        chart1 = CSM.chart(parent, t)
        for (injct, child) in zip(CSM.injections(parent, t), CSM.children(parent, t))
            
            # injct1 = CSM.injections(parent, 1)[1]
            # child1 = chart(fine, CSM.children(parent, 1)[1])

            uv = cartesian(CSM.center(injct))
            @show uv
            p1 = neighborhood(chart1, uv)
            # p2 = CSM.center(child)
            p3 = CSM.center(chart(fine, child))

            # @show x1 = cartesian(p1)
            # @show x2 = cartesian(p2)
            # @show x3 = cartesian(p3)

            # @test isapprox(x1, x2, atol=1e-12)
            @test isapprox(x1, x3, atol=1e-12)
        end
    end
end