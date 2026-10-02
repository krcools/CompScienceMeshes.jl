using Test
using CompScienceMeshes

h = 2π / 51

for T in [Float32, Float64]
    Γ = meshcircle(T(1.0), T(h),2) ## creating the mesh first
    points=vertices(Γ)     ## get the vertices
    num_points=length(points)
    segments=cellarray(Γ)  ## get the faces
    Γ2 = barycentric_refinement(Γ) ## create the refinment of that one
    points2=vertices(Γ2)     ## get the vertices
    segments2=cellarray(Γ2)  ## get the faces

    @test points[1] == points2[1]               # test if the first point in the original and refined mesh are the same
    @test points[num_points]==points2[num_points]# test if the last point in the original has the same index in the refined one
    @test segments[1,1] ==segments2[1,1]        # test if the face are starting with the same point in test
    @test segments[1,2] ==segments2[2,2]        # test if the original second point in the first cell got separated by a single point

    # test barycentric refinement of surface mesh
    local m = meshrectangle(T(1.0), T(1.0), T(0.25), 3)
    local f = barycentric_refinement(m)


    @test CompScienceMeshes.refines(f,m)
    @test numcells(f) == 6*numcells(m)

    local m1 = skeleton(m,1)
    local f1 = skeleton(f,1)
    @test numcells(f1) == 2*numcells(m1) + 6*numcells(m)

    local m0 = skeleton(m,0)
    @test numvertices(f) == numcells(m0) + numcells(m1) + numcells(m)

    ## test bisecting referinment of surfacic meshes
    local b = bisecting_refinement(m)
    @test numcells(b) == 4*numcells(m)

    local b1 = skeleton(b,1)
    @test numcells(b1) == 2*numcells(m1) + 3*numcells(m)
    @test numvertices(b) == numcells(m0) + numcells(m1)

    # the refinement of a curve is a BarycentricRefinement, as for surfaces
    @test CompScienceMeshes.refines(Γ2, Γ)
    @test numcells(Γ2) == 2*numcells(Γ)
    @test numvertices(Γ2) == numvertices(Γ) + numcells(Γ)

    # parent and children are inverse of each other
    local pm = CompScienceMeshes.parent(Γ2)
    for E in 1:numcells(Γ)
        local kids = CompScienceMeshes.children(pm, E)
        @test kids == [2*(E-1)+1, 2*(E-1)+2]
        @test all(CompScienceMeshes.parent(Γ2, k) == E for k in kids)
    end
    # curves in a three-dimensional universe, and with vertex numbering
    # inherited from a parent surface
    for m in (meshsegment(T(1.0), T(1)/4, 3), boundary(meshrectangle(T(1.0), T(1.0), T(1)/2, 3)))
        local f = barycentric_refinement(m)
        @test CompScienceMeshes.refines(f, m)
        @test numcells(f) == 2*numcells(m)
        @test numvertices(f) == numvertices(m) + numcells(m)
    end
    # a submesh of a curve refines too, which the Mesh-only signature could not do
    local Λ = meshsegment(T(1.0), T(1)/4, 3)
    local s = submesh((m,p) -> cartesian(CompScienceMeshes.center(chart(m,p)))[1] < T(1)/2, Λ)
    local sf = barycentric_refinement(s)
    @test CompScienceMeshes.refines(sf, s)
    @test numcells(sf) == 2*numcells(s)
end
