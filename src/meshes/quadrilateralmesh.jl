

struct QuadMesh{T} <: AbstractMesh{3,4,T}
    vertices::Vector{SVector{3,T}}
    faces::Vector{SVector{4,Int}}
end

function indices(m::QuadMesh, cellptr) m.faces[cellptr] end
function vertextype(m::QuadMesh) eltype(m.vertices) end
function dimension(m::QuadMesh) 2 end

function cells(m::QuadMesh) m.faces end
function vertices(m::QuadMesh) m.vertices end

function Base.length(m::QuadMesh) length(m.faces) end
function Base.iterate(m::QuadMesh, state=0) iterate(eachindex(m.faces), state) end

function indextype(m::QuadMesh) SVector{4,Int} end
function indextype(m::QuadMesh, ::Type{Val{N}}) where {N} SVector{N+1,Int} end
function celltype(m::QuadMesh, ::Type{Val{0}}) SimplexGraph{1} end
function celltype(m::QuadMesh, ::Type{Val{1}}) SimplexGraph{2} end
function celltype(m::QuadMesh, ::Type{Val{2}}) QuadrilateralGraph end
function celltype(m::QuadMesh) QuadrilateralGraph end

function chart(m::QuadMesh, cellptr)
    verts = m.vertices[m.faces[cellptr]]
    Quadrilateral(verts...)
end

@testitem "QuadMesh" begin
    using StaticArrays
    p1 = point(0,0,0)
    p2 = point(1,0,0)
    p3 = point(1,1,0)
    p4 = point(0,1,0)
    i1 = SVector(1,2,3,4)
    m = CompScienceMeshes.QuadMesh([p1,p2,p3,p4], [i1])
    @test CompScienceMeshes.indices(m,1) == i1
    @test vertextype(m) == typeof(p1)
    @test dimension(m) == 2
    @test length(m) == 1
    @test CompScienceMeshes.indices(m, first(m)) == i1
end

@testitem "Quadmesh: charts" begin
    using StaticArrays
    p1 = point(0,0,0)
    p2 = point(1,0,0)
    p3 = point(1,1,0)
    p4 = point(0,1,0)
    i1 = SVector(1,2,3,4)
    m = CompScienceMeshes.QuadMesh([p1,p2,p3,p4], [i1])
    ch = chart(m, first(m))
    p = neighborhood(ch, point(0.5, 0.5))
    @test cartesian(p) ≈ point(0.5, 0.5, 0.0) 
end

@testitem "QuadMesh: skeleton" begin
    m = CompScienceMeshes.meshrectangle(2.0, 2.0, 1.0; structured=:quadrilateral)
    @test length(m) == 4

    e = CompScienceMeshes.skeleton_fast(m,Val{1})
    @test length(e) == 12
end

@testitem "QuadMesh: boundary" begin
    m = CompScienceMeshes.meshrectangle(2.0, 2.0, 1.0; structured=:quadrilateral)
    b = CompScienceMeshes.boundary(m)
    @test length(b) == 8
end

@testitem "QuadMesh: connectivity" begin
    m = CompScienceMeshes.meshrectangle(2.0, 2.0, 1.0; structured=:quadrilateral)
    e = CompScienceMeshes.skeleton(m, 1)
    c = CompScienceMeshes.connectivity(e, m, identity)
    @test size(c) == (4,12)
end


function vertexarray(m::QuadMesh)
    [v[i] for v in m.vertices, i in 1:length(eltype(m.vertices))]
end
function cellarray(m::QuadMesh)
    [k[i] for k in m.faces, i in 1:4]
end