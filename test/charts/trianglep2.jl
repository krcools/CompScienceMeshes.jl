@testitem "Quadratic charts: parameter convention" begin

    const CSM = CompScienceMeshes

    p1 = point(1,0,0)
    p2 = point(0,1,0)
    p3 = point(0,0,0)
    p4 = point(0.5,0,0)
    p5 = point(0.5,0.5,0)
    p6 = point(0,0.5,0)

    ch = CSM.TriangleP2(p1,p2,p3,p4,p5,p6)
    u = (1/3, 1/3)
    p = neighborhood(ch, u)
    x = cartesian(p)

    @test x ≈ point(1/3, 1/3, 0)
 
    @test cartesian(neighborhood(ch, (0,0))) ≈ p3
    @test cartesian(neighborhood(ch, (1,0))) ≈ p1
    @test cartesian(neighborhood(ch, (0,1))) ≈ p2

    @test cartesian(neighborhood(ch, (0.5,0))) ≈ p6
    @test cartesian(neighborhood(ch, (0.5,0.5))) ≈ p4
    @test cartesian(neighborhood(ch, (0,0.5))) ≈ p5

end

@testitem "Quadratic charts: faces" begin

    const CSM = CompScienceMeshes

    p1 = point(1,0,0)
    p2 = point(0,1,0)
    p3 = point(0,0,0)
    p4 = point(0.5,0,0)
    p5 = point(0.5,0.5,0)
    p6 = point(0,0.5,0)

    ch = CSM.TriangleP2(p1,p2,p3,p4,p5,p6)

    fc = CSM.faces(ch)

    q1 = neighborhood(fc[1], (0.5,))
    q2 = neighborhood(fc[2], (0.5,))
    q3 = neighborhood(fc[3], (0.5,))

    @test cartesian(q1) ≈ p6
    @test cartesian(q2) ≈ p4
    @test cartesian(q3) ≈ p5

    @test q1.tangents[1] ≈ p2 - p3
    @test q2.tangents[1] ≈ p3 - p1
    @test q3.tangents[1] ≈ p1 - p2

end