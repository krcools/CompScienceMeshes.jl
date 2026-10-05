using Test
using CompScienceMeshes
using StaticArrays

p = simplex(
    point(0,0,0),
    point(1,0,0),
    point(0,1,0))

q1 = simplex(
    point(0.6, 0.6, 0.0),
    point(1.6, 0.6, 0.0),
    point(0.6, 1.6, 0.0),
)

q2 = simplex(
    point(0.4, 0.4, 0.0),
    point(1.4, 0.4, 0.0),
    point(0.4, 1.4, 0.0),
)

@test overlap(p, q1) == false
@test overlap(p, q2) == true

q1 = simplex(
    point(0.949726718617726,-0.278864925637586,0.142314838273285),
    point(0.989821441880933,0.0,-0.142314838273285),
    point(0.989821441880933,0.0,0.142314838273285),
)

q2 = simplex(
    point(0.949726718617726,-0.278864925637586,-0.142314838273285),
    point(0.989821441880933,0.0,-0.142314838273285),
    point(0.949726718617726,-0.278864925637586,0.142314838273285),
)

@test overlap(q1, q2) == false

## make sure the submesh function work for 1D meshes


l1 = meshsegment(1.0,1/2)
vt = skeleton(l1,0)
bd = boundary(l1)

overlaps = overlap_gpredicate(bd)
# pred1 = c -> overlaps(simplex(vertices(vt,c)))
pred1 = c -> overlaps(chart(vt,c))
@test pred1(CompScienceMeshes.SimplexGraph(1))
@test !pred1(CompScienceMeshes.SimplexGraph(2))
@test pred1(CompScienceMeshes.SimplexGraph(3))


# test a case where the segments are:
#   not of unit length
#   colinear and opposite
#   meet in a common point
ch1 = simplex(point(1/3,0,0), point(1/3,1/3,0))
ch2 = simplex(point(1/3,1/3,0), point(1/3,2/3,0))
@test !overlap(ch1, ch2)

@testset "overlap symmetry and boundary cases" begin
    a = simplex(point(0.0, 0.0, 0.0), point(1.0, 0.0, 0.0), point(0.0, 1.0, 0.0))
    identical = a
    contained = simplex(
        point(0.1, 0.1, 0.0), point(0.2, 0.1, 0.0), point(0.1, 0.2, 0.0))
    disjoint = simplex(
        point(2.0, 2.0, 0.0), point(3.0, 2.0, 0.0), point(2.0, 3.0, 0.0))
    touching = simplex(
        point(1.0, 0.0, 0.0), point(2.0, 0.0, 0.0), point(1.0, 1.0, 0.0))

    for q in (identical, contained, disjoint, touching)
        @test overlap(a, q) == overlap(q, a)
    end
    @test overlap(a, identical)
    @test overlap(a, contained)
    @test !overlap(a, disjoint)

    segment = simplex(point(0.0, 0.0, 0.0), point(1.0, 0.0, 0.0))
    segment_overlap = simplex(point(0.5, 0.0, 0.0), point(1.5, 0.0, 0.0))
    segment_touch = simplex(point(1.0, 0.0, 0.0), point(2.0, 0.0, 0.0))
    segment_disjoint = simplex(point(2.0, 0.0, 0.0), point(3.0, 0.0, 0.0))
    @test overlap(segment, segment_overlap)
    @test !overlap(segment, segment_touch)
    @test !overlap(segment, segment_disjoint)
    @test overlap(segment, segment_overlap) == overlap(segment_overlap, segment)
end

@testset "overlap across geometric scales" begin
    for scale in (1e-6, 1.0, 1e6)
        separation = max(10scale, 10sqrt(eps(Float64)))
        p = simplex(
            point(scale, 0.0, 0.0),
            point(0.0, scale, 0.0),
            point(0.0, 0.0, scale),
        )
        q = simplex(
            point(scale, 0.0, 0.0),
            point(0.0, scale, 0.0),
            point(0.0, 0.0, scale),
        )
        far = simplex(
            point(scale + separation, 0.0, 0.0),
            point(0.0, scale + separation, 0.0),
            point(0.0, 0.0, scale + separation),
        )
        @test overlap(p, q)
        @test !overlap(p, far)
        @test overlap(p, far) == overlap(far, p)
    end
end

@testset "point overlap tolerance scales with magnitude" begin
    for scale in (1.0, 1e6)
        shift = scale
        p = simplex(point(shift, shift, shift))
        q = simplex(point(shift, shift, shift))
        inside = simplex(point(shift + 0.5sqrt(eps(Float64)) * scale, shift, shift))
        outside = simplex(point(shift + 2sqrt(eps(Float64)) * scale, shift, shift))
        @test overlap(p, q)
        @test overlap(p, inside)
        @test !overlap(p, outside)
        @test overlap(p, outside) == overlap(outside, p)
    end
end

@testset "small positive overlap across scales" begin
    for scale in (1e-6, 1.0, 1e6)
        delta = 1e-3
        shift = (1 - delta) * scale

        segment = simplex(point(0.0, 0.0, 0.0), point(scale, 0.0, 0.0))
        shifted_segment = simplex(point(shift, 0.0, 0.0), point(shift + scale, 0.0, 0.0))
        @test overlap(segment, shifted_segment)

        triangle = simplex(
            point(0.0, 0.0, 0.0), point(scale, 0.0, 0.0), point(0.0, scale, 0.0))
        shifted_triangle = simplex(
            point(shift, 0.0, 0.0), point(shift + scale, 0.0, 0.0),
            point(shift, scale, 0.0))
        @test overlap(triangle, shifted_triangle)

        tetrahedron = simplex(
            point(0.0, 0.0, 0.0), point(scale, 0.0, 0.0),
            point(0.0, scale, 0.0), point(0.0, 0.0, scale))
        # A vertex of the shifted tetrahedron lies a small distance inside the
        # first one (vertex-in-interior is the tetrahedron overlap criterion).
        offset = delta * scale
        shifted_tetrahedron = simplex(
            point(offset, offset, offset), point(offset + scale, offset, offset),
            point(offset, offset + scale, offset), point(offset, offset, offset + scale))
        @test overlap(shifted_tetrahedron, tetrahedron)
        @test overlap(tetrahedron, shifted_tetrahedron)
    end
end

@testset "tetrahedron overlap across scales and translations" begin
    for scale in (1e-6, 1.0, 1e6)
        shift = point(3scale, -2scale, scale)
        p = simplex(
            shift + point(0.0, 0.0, 0.0),
            shift + point(scale, 0.0, 0.0),
            shift + point(0.0, scale, 0.0),
            shift + point(0.0, 0.0, scale),
        )
        q = simplex(
            shift + point(0.1scale, 0.1scale, 0.1scale),
            shift + point(0.2scale, 0.1scale, 0.1scale),
            shift + point(0.1scale, 0.2scale, 0.1scale),
            shift + point(0.1scale, 0.1scale, 0.2scale),
        )
        # q lies inside p: a vertex of q is in the interior of p.
        @test overlap(q, p)
    end
end

@testset "tetrahedron overlap criterion and boundary semantics" begin
    # overlap(p, q) is true when a vertex of p, or the center of p, lies in the
    # interior of q. It is therefore not symmetric, and contact does not count.
    p = simplex(
        point(0.0, 0.0, 0.0),
        point(2.0, 0.0, 0.0),
        point(0.0, 2.0, 0.0),
        point(0.0, 0.0, 2.0),
    )
    inside = simplex(
        point(0.1, 0.1, 0.1),
        point(0.5, 0.1, 0.1),
        point(0.1, 0.5, 0.1),
        point(0.1, 0.1, 0.5),
    )
    @test overlap(inside, p)
    @test overlap(p, p)

    shared_vertex = simplex(
        point(2.0, 0.0, 0.0),
        point(3.0, 0.0, 0.0),
        point(2.0, 1.0, 0.0),
        point(2.0, 0.0, 1.0),
    )
    @test !overlap(p, shared_vertex)
    @test !overlap(shared_vertex, p)

    shared_face = simplex(
        point(0.0, 0.0, 0.0),
        point(0.0, 2.0, 0.0),
        point(0.0, 0.0, 2.0),
        point(-1.0, 0.0, 0.0),
    )
    @test !overlap(p, shared_face)
    @test !overlap(shared_face, p)

    separated = simplex(
        point(2.5, 0.0, 0.0),
        point(3.5, 0.0, 0.0),
        point(2.5, 1.0, 0.0),
        point(2.5, 0.0, 1.0),
    )
    @test !overlap(p, separated)
    @test !overlap(separated, p)
end

@testset "triangle overlap symmetry and invariance" begin
    for scale in (1e-6, 1.0, 1e6)
        shift = point(2scale, -scale, 3scale)
        p = simplex(
            shift + point(0.0, 0.0, 0.0),
            shift + point(scale, 0.0, 0.0),
            shift + point(0.0, scale, 0.0),
        )
        contained = simplex(
            shift + point(0.1scale, 0.1scale, 0.0),
            shift + point(0.2scale, 0.1scale, 0.0),
            shift + point(0.1scale, 0.2scale, 0.0),
        )
        separated = simplex(
            shift + point(2scale, 2scale, 0.0),
            shift + point(3scale, 2scale, 0.0),
            shift + point(2scale, 3scale, 0.0),
        )

        @test overlap(p, contained)
        @test overlap(contained, p)
        @test overlap(p, contained) == overlap(contained, p)
        @test !overlap(p, separated)
        @test !overlap(separated, p)
    end

    scale = 1e6
    depth = 1e-6 * scale
    p = simplex(
        point(0.0, 0.0, 0.0),
        point(scale, 0.0, 0.0),
        point(0.0, scale, 0.0),
    )
    slightly_overlapping = simplex(
        point(scale - depth, 0.0, 0.0),
        point(2scale - depth, 0.0, 0.0),
        point(scale - depth, scale, 0.0),
    )
    @test overlap(p, slightly_overlapping)
    @test overlap(slightly_overlapping, p)
end

@testset "degenerate overlap inputs" begin
    point_a = simplex(point(1.0, 2.0, 3.0))
    point_b = simplex(point(1.0, 2.0, 3.0))
    point_c = simplex(point(1.0 + 1e-6, 2.0, 3.0))
    @test overlap(point_a, point_b)
    @test !overlap(point_a, point_c)

    zero_segment = simplex(point(0.0, 0.0, 0.0), point(0.0, 0.0, 0.0))
    segment = simplex(point(0.0, 0.0, 0.0), point(1.0, 0.0, 0.0))
    @test !overlap(zero_segment, segment)
    @test !overlap(segment, zero_segment)
end

@testset "overlap predicates do not allocate" begin
    point_patch = simplex(point(0.0, 0.0, 0.0))
    segment_patch = simplex(point(0.0, 0.0, 0.0), point(1.0, 0.0, 0.0))
    triangle_patch = simplex(
        point(0.0, 0.0, 0.0), point(1.0, 0.0, 0.0), point(0.0, 1.0, 0.0))
    tetrahedron_patch = simplex(
        point(0.0, 0.0, 0.0), point(1.0, 0.0, 0.0),
        point(0.0, 1.0, 0.0), point(0.0, 0.0, 1.0))

    overlap(point_patch, point_patch)
    overlap(segment_patch, segment_patch)
    overlap(triangle_patch, triangle_patch)
    overlap(tetrahedron_patch, tetrahedron_patch)

    @test @allocated(overlap(point_patch, point_patch)) == 0
    @test @allocated(overlap(segment_patch, segment_patch)) == 0
    @test @allocated(overlap(triangle_patch, triangle_patch)) == 0
    @test @allocated(overlap(tetrahedron_patch, tetrahedron_patch)) == 0
end

@testset "realistic tetrahedral mesh overlap" begin
    mesh = CompScienceMeshes.read_gmsh3d_mesh(
        joinpath(@__DIR__, "assets", "cuboid_tetrahedron.msh"))
    sample_ids = 1:max(1, cld(numcells(mesh), 100)):numcells(mesh)

    @test all(i -> begin
        cell = chart(mesh, i)
        overlap(cell, cell)
    end, sample_ids)

    # Neighboring cells of a conforming mesh only touch, which is not overlap.
    @test all(i -> begin
        j = mod1(i + 1, numcells(mesh))
        p = chart(mesh, i)
        q = chart(mesh, j)
        !overlap(p, q) && !overlap(q, p)
    end, sample_ids)

    p = chart(mesh, first(sample_ids))
    c = 0.7 * p.vertices[1] + 0.1 * p.vertices[2] +
        0.1 * p.vertices[3] + 0.1 * p.vertices[4]
    contained = simplex((c + 0.1 * (v - c) for v in p.vertices)...)
    @test overlap(contained, p)

    far = simplex(
        p.vertices[1] + point(100.0, 0.0, 0.0),
        p.vertices[2] + point(100.0, 0.0, 0.0),
        p.vertices[3] + point(100.0, 0.0, 0.0),
        p.vertices[4] + point(100.0, 0.0, 0.0),
    )
    @test !overlap(p, far)
    @test !overlap(far, p)
end

@testset "large tetrahedral star overlap" begin
    mesh = CompScienceMeshes.readmesh(
        joinpath(@__DIR__, "assets", "tetra_star.in"))
    @test numcells(mesh) == 12288

    @test all(i -> begin
        cell = chart(mesh, i)
        overlap(cell, cell)
    end, 1:numcells(mesh))

    sample_ids = 1:max(1, cld(numcells(mesh), 1000)):numcells(mesh)
    @test all(i -> begin
        j = mod1(i + 1, numcells(mesh))
        p = chart(mesh, i)
        q = chart(mesh, j)
        !overlap(p, q) && !overlap(q, p)
    end, sample_ids)

    p = chart(mesh, first(sample_ids))
    far = simplex(
        p.vertices[1] + point(100.0, 0.0, 0.0),
        p.vertices[2] + point(100.0, 0.0, 0.0),
        p.vertices[3] + point(100.0, 0.0, 0.0),
        p.vertices[4] + point(100.0, 0.0, 0.0),
    )
    @test !overlap(p, far)
    @test !overlap(far, p)
end

@testset "large-scale segment mesh" begin
    sensitive = simplex(point(-0.5, 0.0, 0.0), point(0.5, 0.0, 0.0))
    short = simplex(point(-0.5, 0.0, 0.0), point(-0.499, 1e-13, 0.0))
    @test overlap(sensitive, short)
    @test overlap(short, sensitive)

    # A star of 100000 segments radiating from the origin, generated here
    # rather than read from a large fixture file.
    nsegments = 100000
    vertices = Vector{typeof(point(0.0, 0.0, 0.0))}(undef, 2 * nsegments)
    faces = Vector{CompScienceMeshes.SimplexGraph{2}}(undef, nsegments)
    for i in 1:nsegments
        theta = 2pi * i / nsegments
        phi = pi * mod(i * 0.618033988749895, 1.0)
        direction = point(sin(phi) * cos(theta), sin(phi) * sin(theta), cos(phi))
        vertices[2i - 1] = (0.1 + mod(i, 7)) * direction
        vertices[2i] = (1.1 + mod(i, 7)) * direction
        faces[i] = CompScienceMeshes.SimplexGraph(2i - 1, 2i)
    end
    mesh = Mesh(vertices, faces)
    @test numcells(mesh) == nsegments

    @test all(i -> begin
        segment = chart(mesh, i)
        overlap(segment, segment)
    end, 1:numcells(mesh))

    # Segments on different rays only touch (at most) and must not overlap.
    @test all(i -> begin
        !overlap(chart(mesh, i), chart(mesh, mod1(i + 1, numcells(mesh))))
    end, 1:100:numcells(mesh))
end

@testset "realistic sphere overlap and predicates" begin
    sphere = readmesh(joinpath(@__DIR__, "assets", "sphere8.in"))
    overlap_predicate = overlap_gpredicate(sphere)
    inclosure_predicate = inclosure_gpredicate(sphere)

    for cell in cells(sphere)
        patch = chart(sphere, cell)
        @test overlap(patch, patch)
        @test overlap_predicate(patch)
    end

    for cell in cells(sphere), vertex in vertices(chart(sphere, cell))
        @test inclosure_predicate(vertex)
    end
    @test !inclosure_predicate(point(0.0, 0.0, 0.0))
    @test !inclosure_predicate(point(10.0, 0.0, 0.0))

    cells_list = collect(cells(sphere))
    neighboring_pair = nothing
    for i in eachindex(cells_list), j in (i+1):length(cells_list)
        p = chart(sphere, cells_list[i])
        q = chart(sphere, cells_list[j])
        shared_vertices = count(v -> any(v == w for w in q.vertices), p.vertices)
        if shared_vertices >= 2
            neighboring_pair = (p, q)
            break
        end
    end
    @test neighboring_pair !== nothing
    if neighboring_pair !== nothing
        p, q = neighboring_pair
        @test overlap(p, q) == overlap(q, p)
    end

    first_patch = chart(sphere, first(cells(sphere)))
    far_patch = simplex(
        first_patch.vertices[1] + point(10.0, 0.0, 0.0),
        first_patch.vertices[2] + point(10.0, 0.0, 0.0),
        first_patch.vertices[3] + point(10.0, 0.0, 0.0),
    )
    @test !overlap(first_patch, far_patch)
    @test overlap(first_patch, far_patch) == overlap(far_patch, first_patch)
    @test !overlap_predicate(far_patch)
end
