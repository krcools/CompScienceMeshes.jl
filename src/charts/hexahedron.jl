struct Hexahedron{P}
    vertices::SVector{8, P}
    p21::P   # p2 - p1
    p41::P   # p4 - p1
    p51::P   # p5 - p1
    p1243::P # p1 - p2 - p4 + p3
    p1256::P # p1 - p2 - p5 + p6
    p1458::P # p1 - p4 - p5 + p8 
    p1to8::P #-p1 + p2 - p3 + p4 + p5 - p6 + p7 - p8
end

Hexahedron(p1::P, p2::P, p3::P, p4::P, p5::P, p6::P, p7::P, p8::P) where P = Hexahedron{P}(
    SVector(p1, p2, p3, p4, p5, p6, p7, p8), 
    p2 - p1,
    p4 - p1, 
    p5 - p1,
    p1 - p2 - p4 + p3,
    p1 - p2 - p5 + p6,
    p1 - p4 - p5 + p8,
    -p1 + p2 - p3 + p4 + p5 - p6 + p7 - p8
)

Hexahedron(p) = Hexahedron(p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8])

function coordtype(hex::Hexahedron{P}) where P eltype(P) end

function vertices(hex::Hexahedron) hex.vertices end

function cartesian(hex::Hexahedron, u)
    u1, u2, u3 = u

    return hex.vertices[1] + u1*hex.p21 + u2*hex.p41 + u3*hex.p51 + 
    u1*u2*hex.p1243 +
    u1*u3*hex.p1256 +
    u2*u3*hex.p1458 +
    u1*u2*u3*hex.p1to8
end

function tangents(hex::Hexahedron, u)
    u1, u2, u3 = u

    ∂ru1 = hex.p21 + u2*hex.p1243 + u3*hex.p1256 + u2*u3*hex.p1to8
    ∂ru2 = hex.p41 + u1*hex.p1243 + u3*hex.p1458 + u1*u3*hex.p1to8
    ∂ru3 = hex.p51 + u1*hex.p1256 + u2*hex.p1458 + u1*u2*hex.p1to8

    return SMatrix{3,3}(
        ∂ru1[1], ∂ru1[2], ∂ru1[3], # ---> column 1 for an SMatrix
        ∂ru2[1], ∂ru2[2], ∂ru2[3], # ---> column 2 for an SMatrix
        ∂ru3[1], ∂ru3[2], ∂ru3[3]  # ---> column 3 for an SMatrix
    )
end

function jacobian(hex::Hexahedron, u)

    J = tangents(hex, u)

    return dot(cross(J[:,1], J[:,2]), J[:,3])
end

function jacobian_(hex::Hexahedron, ∂ru1, ∂ru2, ∂ru3)

    return dot(cross(∂ru1, ∂ru2), ∂ru3)
end




struct NeighborhoodHex{C,P,Q,T,J}
    chart::C
    parametric::P
    cartesian::Q
    tangents::T
    jacobian::J
end

function neighborhood(hex::Hexahedron, u)
    c = cartesian(hex, u)
    J = tangents(hex, u)
    j = jacobian_(hex, J[:,1], J[:,2], J[:,3])
    
    return NeighborhoodHex(hex, u, c, J, j)
end

function parametric(mp::NeighborhoodHex) mp.parametric end
function cartesian(mp::NeighborhoodHex) mp.cartesian end
function tangents(mp::NeighborhoodHex) mp.tangents end
function tangents(mp::NeighborhoodHex, i::Int) mp.tangents[:,i] end
function jacobian(mp::NeighborhoodHex) mp.jacobian end




struct RefHexahedron{T} end

domain(hex::Hexahedron{P}) where {P} = RefHexahedron{eltype(P)}()

function vertices(hex::RefHexahedron{T}) where T
    SVector(
        point(T, 0, 0, 0), #p1 , A
        point(T, 1, 0, 0), #p2 , B
        point(T, 1, 1, 0), #p3 , C
        point(T, 0, 1, 0), #p4 , D
        point(T, 0, 0, 1), #p5 , E
        point(T, 1, 0, 1), #p6 , F
        point(T, 1, 1, 1), #p7 , G
        point(T, 0, 1, 1)  #p8 , H
        )
end

function permute_vertices(hex::RefHexahedron, I)
    V = vertices(hex)
    return Hexahedron(V[I[1]], V[I[2]], V[I[3]], V[I[4]], V[I[5]], V[I[6]], V[I[7]], V[I[8]])
end





@testitem "hexahedron" begin

    # cube
    p1 = point(0,0,0)
    p2 = point(1,0,0)
    p3 = point(1,1,0)
    p4 = point(0,1,0)
    p5 = point(0,0,1)
    p6 = point(1,0,1)
    p7 = point(1,1,1)
    p8 = point(0,1,1)

    hex = CompScienceMeshes.Hexahedron(p1,p2,p3,p4,p5,p6,p7,p8)

    @test coordtype(hex) == Float64

    mp = neighborhood(hex, (0.5,0.5,0.5))
    @test cartesian(mp) ≈ point(0.5,0.5,0.5)
    @test tangents(mp, 1) ≈ point(1,0,0)
    @test tangents(mp, 2) ≈ point(0,1,0)
    @test tangents(mp, 3) ≈ point(0,0,1)
    @test hex.p1243 ≈ point(0,0,0)
    @test hex.p1256 ≈ point(0,0,0)
    @test hex.p1458 ≈ point(0,0,0)
    @test hex.p1to8 ≈ point(0,0,0)


    # distorted hexahedron
    p1 = point(0.0, 0.0, 0.0)
    p2 = point(2.2, 0.6, -0.9)
    p3 = point(2.8, 2.1, -0.2)
    p4 = point(-1.0, 0.7, -0.6)
    p5 = point(0.0, 0.0, 0.9)
    p6 = point(2.3, -0.2, 5.2)
    p7 = point(2.7, 1.1, 3.6)
    p8 = point(-0.2, 2.2, 3.1)

    hex = CompScienceMeshes.Hexahedron(p1,p2,p3,p4,p5,p6,p7,p8)

    @test cartesian(neighborhood(hex, (0,0,0))) ≈ p1
    @test cartesian(neighborhood(hex, (1,0,0))) ≈ p2
    @test cartesian(neighborhood(hex, (1,1,0))) ≈ p3 
    @test cartesian(neighborhood(hex, (0,1,0))) ≈ p4 
    @test cartesian(neighborhood(hex, (0,0,1))) ≈ p5 
    @test cartesian(neighborhood(hex, (1,0,1))) ≈ p6
    @test cartesian(neighborhood(hex, (1,1,1))) ≈ p7 
    @test cartesian(neighborhood(hex, (0,1,1))) ≈ p8
    mp = neighborhood(hex, (0,0,0))
    @test tangents(mp, 1) ≈ hex.p21
    @test tangents(mp, 2) ≈ hex.p41
    @test tangents(mp, 3) ≈ hex.p51
    @test cartesian(neighborhood(hex, (0.5,0.5,0.5))) ≈ sum(hex.vertices)/8


    # reference hexahedron
    refhex = domain(hex)
    @test typeof(refhex).parameters[1] == Float64
    @test vertices(refhex)[5] == point(Float64, 0, 0, 1)
end