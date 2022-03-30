using Test
using CompScienceMeshes
using LinearAlgebra
using StaticArrays

for T in [Float64, Float32]
    rectangle = meshrectangle(T(1.0), T(1.0), T(1.0));

    @test vertices(rectangle,1) == point(T, 0.0, 0.0, 0.0)
    @test vertices(rectangle,2) == point(T, 0.0, 1.0, 0.0)
    @test vertices(rectangle,3) == point(T, 1.0, 0.0, 0.0)

    @test numvertices(rectangle) == 4
    @test numcells(rectangle) == 2
    @test dimension(rectangle) == 2

    faces = skeleton(rectangle, 2)
    edges = skeleton(rectangle, 1)
    verts = skeleton(rectangle, 0)

    @test numcells(verts) == 4
    @test numcells(edges) == 5
    @test numcells(faces) == 2

    Λ  = connectivity(verts, edges, sign)
    Σᵀ = connectivity(edges, faces, sign)

    @test size(Λ)  == (5,4)
    @test size(Σᵀ) == (2,5)

    @test norm(Σᵀ*Λ, Inf) ≈ 0

    bnd = boundary(rectangle)

    @test dimension(bnd) == 1
    @test numvertices(bnd) == 4
    @test numcells(bnd) == 4

    vtoe, num_vtoe = vertextocellmap(edges)
    vtof, num_vtof = vertextocellmap(faces)

    @test size(vtoe) == (4,3)
    @test size(vtof) == (4,2)

    @test size(num_vtoe) == (4,)
    @test size(num_vtof) == (4,)

    file = joinpath(dirname(@__FILE__), "sphere.in")
    sphere = readmesh(file)
    @test numvertices(sphere) == 335
    @test numcells(sphere) == 666

    m = meshrectangle(T(1.0),T(1.0),T(1.0))
    CompScienceMeshes.rotate!(m, [1,0,0]*T(π/2))

    @test CompScienceMeshes.isoriented(boundary(m))

    tol = sqrt(eps(T))
    @test norm(vertices(m,1) - point(T,0.0,0.0,0.0)) < tol
    @test norm(vertices(m,2) - point(T,0.0,0.0,1.0)) < tol
    @test norm(vertices(m,3) - point(T,1.0,0.0,0.0)) < tol
    @test norm(vertices(m,4) - point(T,1.0,0.0,1.0)) < tol

    ## test fliporientation
    V = [
        point(T,0,0,0),
        point(T,1,0,0),
        point(T,0,1,0),
    ]

    F = [
        index(1,2,3)
    ]

    m = Mesh(V,F)
    n = fliporientation(m)

    @test m.faces[1] == index(1,2,3)
    @test n.faces[1] == index(2,1,3)

    # test the construction of transposed connectivity matricees
    Σᵀ = connectivity(edges, faces, identity)
    Σ = connectivity(faces, edges, identity)
    @test norm(Σᵀ - Σ', Inf) == 0
end