using Test
using CompScienceMeshes
using StaticArrays
using LinearAlgebra

# Locating points in curvilinear triangles, and the overlap and intersection of a cell of a mesh
# with the cells of its barycentric refinement: exact for a cell and its children, empty for the
# children of its neighbours, which share at most an edge with it.
@testset "intersection of curvilinear triangles, order $O" for O in 2:3

    CSM = CompScienceMeshes
    Γ = gmshsphere(1.0, 0.6; order=O)
    R = barycentric_refinement(Γ)

    for c in 1:numcells(Γ)
        P = chart(Γ, c)

        # carttobary inverts barytocart, also on the continuation beyond the reference triangle
        for u in (SVector(0.2, 0.3), SVector(0.7, 0.1), SVector(-0.3, 0.6), SVector(1.2, 0.4))
            @test carttobary(P, barytocart(P, u)) ≈ u atol=1e-12
        end

        for k in 6(c-1)+1:6c
            Q = chart(R, k)
            @test overlap(P, Q)
            @test overlap(Q, P)
            for (A, B) in ((P, Q), (Q, P))
                pieces = intersection(A, B)
                @test sum(volume, pieces) ≈ volume(Q) rtol=1e-12
                # the pieces have the orientation of the first operand
                @test all(dot(normal(center(s)), normal(center(A))) > 0 for s in pieces)
            end
            for kk in 6(c-1)+1:6c
                kk == k || @test !overlap(Q, chart(R, kk))
            end
        end

        for d in 1:numcells(Γ)
            d == c && continue
            isempty(intersect(CSM.indices(Γ, c), CSM.indices(Γ, d))) && continue
            for k in 6(d-1)+1:6d
                @test !overlap(P, chart(R, k))
                @test isempty(intersection(P, chart(R, k)))
            end
        end
    end
end
