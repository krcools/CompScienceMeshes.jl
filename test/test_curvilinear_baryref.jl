using Test
using CompScienceMeshes
using StaticArrays
using LinearAlgebra
using SparseArrays

const CSM = CompScienceMeshes

@testset "barycentric refinement of curvilinear triangles, order $O" for O in 1:3

    # gmsh returns a flat Mesh for order 1: the curvilinear mesh of order 1 is built from it
    Γ = if O == 1
        Γflat = gmshsphere(1.0, 0.6)
        CSM.CurvilinearMesh(vertices(Γflat), [CSM.HOSimplexGraph{3,1}(c.indices...) for c in cells(Γflat)])
    else
        gmshsphere(1.0, 0.6; order=O)
    end
    @test Γ isa CSM.CurvilinearMesh
    @test CSM.indextype(Γ) == SVector{3,Int}

    # the center of a curved triangle is the image of the reference centroid
    for p in 1:numcells(Γ)
        ch = chart(Γ, p)
        @test parametric(center(ch)) ≈ SVector(1 / 3, 1 / 3)
        @test cartesian(center(ch)) ≈ barytocart(ch, SVector(1 / 3, 1 / 3))
    end
    O == 1 && @test all(cartesian(center(chart(Γ, p))) ≈ cartesian(center(chart(Γflat, p))) for p in 1:numcells(Γ))

    R = barycentric_refinement(Γ)
    F = R.mesh

    @test R isa CSM.BarycentricRefinement
    @test CSM.refines(R, Γ)
    @test F isa CSM.CurvilinearMesh
    @test CSM.meshorder(F) == O
    @test numcells(R) == 6 * numcells(Γ)

    # the children reproduce their parent exactly: a child's reference point μ is the parent's
    # point Σ μ_i Λ_i, with Λ_i the parent-barycentric coordinates of the child's corners
    corners = (SVector(1.0, 0, 0), SVector(0, 1.0, 0), SVector(0, 0, 1.0))
    ctr = SVector(1 / 3, 1 / 3, 1 / 3)
    μs = [SVector(0.2, 0.3), SVector(0.6, 0.1), SVector(0.05, 0.9), SVector(1 / 3, 1 / 3)]
    maxdev = 0.0
    for p in 1:numcells(Γ)
        chp = chart(Γ, p)
        for j in 1:3
            a, b = corners[mod1(j + 1, 3)], corners[mod1(j + 2, 3)]
            m = (a + b) / 2
            for (s, Λ) in enumerate(((a, m, ctr), (b, ctr, m)))
                chc = chart(F, 6 * (p - 1) + 2 * (j - 1) + s)
                for μ in μs
                    λ = μ[1] * Λ[1] + μ[2] * Λ[2] + (1 - μ[1] - μ[2]) * Λ[3]
                    maxdev = max(maxdev, norm(barytocart(chc, μ) - barytocart(chp, SVector(λ[1], λ[2]))))
                end
            end
        end
    end
    @test maxdev < 1e-13

    # conforming: every edge of the refined (closed) surface is shared by exactly two cells, with
    # opposite orientations, and no node appears twice
    edges = skeleton(F, 1)
    D = connectivity(edges, F)
    @test all(i -> length(nzrange(D, i)) == 2, 1:size(D, 2))
    @test all(i -> sum(nonzeros(D)[nzrange(D, i)] .|> sign) == 0, 1:size(D, 2))
    used = unique(reduce(vcat, [collect(c.indices) for c in cells(F)]))
    @test length(unique(map(v -> round.(v, digits=10), vertices(F)[used]))) == length(used)

    # the Euler characteristic of the sphere, from the corners of the refined cells
    nv = length(unique(reduce(vcat, [collect(CSM.indices(c)) for c in cells(F)])))
    @test nv - numcells(edges) + numcells(F) == 2

    # the children keep the orientation of their parent: their normals agree at a common point
    for p in 1:numcells(Γ), k in 1:6
        np = normal(neighborhood(chart(Γ, p), SVector(1 / 3, 1 / 3)))
        chc = chart(F, 6 * (p - 1) + k)
        nc = normal(neighborhood(chc, SVector(1 / 3, 1 / 3)))
        @test dot(np, nc) > 0
    end

    # a flat parent (order 1): the refinement coincides with the one of the flat mesh
    if O == 1
        Rflat = barycentric_refinement(Γflat)
        @test numcells(Rflat) == numcells(R)
        for q in 1:numcells(R)
            chf = chart(Rflat, q)
            chc = chart(R, q)
            @test sort(collect(CSM.vertices(chf)), by=Tuple) ≈ sort(collect(CSM.nodes(chc)), by=Tuple)
        end
    end
end
