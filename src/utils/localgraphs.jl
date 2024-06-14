abstract type AbstractGraph end
function Base.getindex(g::AbstractGraph, i::Int) g.indices[i] end
function Base.hash(g::AbstractGraph, a::UInt) hash(sort(g.indices), a) end
function Base.length(g::AbstractGraph) length(g.indices) end

struct QuadrilateralGraph <: AbstractGraph
    indices::SVector{4,Int}
end

function dimension(g::QuadrilateralGraph) 2 end
function dimtype(g::QuadrilateralGraph) Val{2} end

struct SimplexGraph{N} <: AbstractGraph
    indices::SVector{N,Int}
end

function dimension(h::SimplexGraph{N}) where {N} N-1 end
function dimtype(g::SimplexGraph{N}) where {N} Val{N-1} end

function Base.isequal(g::G, h::G) where {G <: AbstractGraph}
    p = indexin(g.indices, h.indices)
    any(p .== nothing) && return false
    return true
end

function skeleton(g::SimplexGraph{N}, ::Type{Val{M}}) where {N,M}
    collect(idcs for idcs in combinations(g.indices, M+1))
end

# function skeleton(g::SimplexGraph{N}, ::Type{Val{N-1}}) where {N}
#     SVector(g.indices)
# end

function skeleton(g::SimplexGraph{3}, ::Type{Val{1}})
    return SVector(
        (g[2],g[3]),
        (g[3],g[1]),
        (g[1],g[2]),
    )
end


@testitem "skeleton SimplexGraph" begin
    g = CompScienceMeshes.SimplexGraph{4}((1,2,3,4))
    edges = skeleton(g,Val{1})
    @test length(edges) == 6
    # @test eltype(edges) == CompScienceMeshes.SVector{2,Int}
end

function skeleton(g::QuadrilateralGraph, ::Type{Val{1}})
    return SVector(
        (g[1],g[2]),
        (g[2],g[3]),
        (g[3],g[4]),
        (g[4],g[1]),)
end

function relorientation(g::G, h::AbstractGraph) where {G <: SimplexGraph}
    dimension(g) > dimension(h) && return relorientation(h, g)
    faces = collect(skeleton(h, dimtype(g)))
    i = something(findfirst(face -> isequal(G(face), g), faces), 0)
    i == 0 && return 0
    p = Permutations.Permutation(indexin(g.indices, SVector(faces[i])))
    σ = Permutations.sign(p)
    return σ * i
end

function relorientation(g::G, h::H) where {G <: SimplexGraph, H <: SimplexGraph}
    dimension(g) > dimension(h) && return relorientation(h, g)
    return relorientation(g.indices, h.indices)
end

# function relorientation(g::G, h::H) where {G <: AbstractGraph, H <: AbstractGraph}
#     dimension(g) > dimension(h) && return relorientation(h, g)
#     faces = collect(skeleton(h, dimtype(g)))
#     i = findfirst(face -> isequal(G(face), g), faces)
#     @assert i != nothing
#     p = Permutations.Permutation(indexin(g.indices, SVector(faces[i])))
#     σ = Permutations.sign(p)
#     return σ * i
# end

@testitem "relorientation Simplex 2D" begin
    h = CompScienceMeshes.SimplexGraph{3}([1,2,3])
    g = CompScienceMeshes.SimplexGraph{2}([3,2])
    σi = CompScienceMeshes.relorientation(g, h)
    @test σi == -1
end

@testitem "relorientation Quadrilateral" begin
    h = CompScienceMeshes.QuadrilateralGraph([1,2,3,4])
    g = CompScienceMeshes.SimplexGraph{2}([3,2])
    σi = CompScienceMeshes.relorientation(g, h)
    @test σi == -2
end