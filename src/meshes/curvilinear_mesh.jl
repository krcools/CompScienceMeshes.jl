# Mesh: 1D elements embedded in R^U, with N control points per element (N = order+1)
mutable struct CurvilinearMesh{U,D1,N,T,O} <: AbstractMesh{U,N,T}
    vertices::Vector{SVector{U,T}}   # coordinates in R^U
    # faces::Vector{SVector{N,Int}}    # connectivity (local order = Gmsh order)
    faces::Vector{HOSimplexGraph{D1,O,N}} # connectivity (local order = Gmsh order)

   
end

function CurvilinearMesh(
    vertices::AbstractVector{SVector{U,T}},
    # faces::AbstractVector{SVector{N,Int}},
    faces::AbstractVector{HOSimplexGraph{D1,O,N}},
    manifolddimension,
    order::Integer
) where {U,N,T,D1,O}
    order > 0 || throw(ArgumentError("order must be positive."))
    order == O || throw(ArgumentError("order inconsistent with cell type."))
    manifolddimension + 1 == D1 || throw(ArgumentError("manifold dimension inconsistent with cell type."))
    CurvilinearMesh{U,D1,N,T,O}(vertices, faces)
end

# This constructor acts like `similar` for arrays. It allows to override the type
# parameters of the mesh, without referring to the abstract parent type.
function similar_mesh(
    mesh::CurvilinearMesh{U,D1,N,T,O},
    vertices::AbstractVector{SVector{U,T}},
    faces::AbstractVector{HOSimplexGraph{D1new,O,Nnew}},
) where {U,N,T,D1,O,D1new,Nnew}
    CurvilinearMesh{U,D1new,Nnew,T,O}(vertices, faces)
end


function CurvilinearMesh(
    vertices::AbstractVector{SVector{U,T}},
    faces::AbstractVector{HOSimplexGraph{D1,O,N}}
) where {U,N,T,D1,O}
    CurvilinearMesh(vertices, faces, D1-1, O)
end




                    #=

# Constructors from abstract containers
CurvilinearMesh(vertices::AbstractVector{<:SVector{U,T}},
                faces::AbstractVector{<:SVector{N,Int}},
                manifolddimension::Integer,
                order::Integer) where {U,N,T} =
    CurvilinearMesh(Vector{SVector{U,T}}(vertices),
                    Vector{SVector{N,Int}}(faces),
                    manifolddimension,
                    order)

# From matrices (verts: U×nV, faces: N×nF)
function CurvilinearMesh(verts::AbstractMatrix{T},
                         faces::AbstractMatrix{<:Integer},
                         manifolddimension,
                         order::Integer) where {T}
    U = size(verts, 1)
    N = size(faces, 1)
    v = [SVector{U,T}(verts[:, i]) for i in 1:size(verts, 2)]
    f = [SVector{N,Int}(faces[:, i]) for i in 1:size(faces, 2)]
    CurvilinearMesh(v, f, order)
end
=#

# Mesh interface
coordtype(::CurvilinearMesh{U,D1,N,T,O}) where {U,D1,N,T,O} = T
vertextype(::CurvilinearMesh{U,D1,N,T,O}) where {U,D1,N,T,O} = SVector{U,T}
universedimension(::CurvilinearMesh{U}) where {U} = U
dimension(::CurvilinearMesh{U,D1}) where {U,D1} = D1-1

function celltype(m::CurvilinearMesh{U,D1,N,T,O}) where {U,D1,N,T,O}
    # SimplexGraph{N-O+1}
    HOSimplexGraph{D1,O,N}
end

function celltype(m::CurvilinearMesh{U,D1,N,T,O}, ::Type{Val{M}}) where {U,D1,N,T,O,M}
    # SimplexGraph{M+1}
    HOSimplexGraph{M+1,O,binomial(M+O,O)}
end

function indextype(m::CurvilinearMesh{U,D1,N,T,O}) where {U,D1,N,T,O} SVector{N-O+1,Int} end
function indextype(m::CurvilinearMesh{U,D1,N,T,O}, ::Type{Val{M}}) where {U,D1,N,T,O,M} SVector{M+1,Int} end

function indices(m::CurvilinearMesh{U,1,N,T,O}, cellptr) where {U,N,T,O}
    cell = m.faces[cellptr]
    return SVector(cell.indices[begin])
end

function indices(m::CurvilinearMesh{U,2,N,T,O}, cellptr) where {U,N,T,O}
    cell = m.faces[cellptr]
    return SVector(cell.indices[begin], cell.indices[begin+1])
end

function indices(m::CurvilinearMesh{U,3,N,T,O}, cellptr) where {U,N,T,O}
    cell = m.faces[cellptr]
    return SVector(cell.indices[begin], cell.indices[begin+1], cell.indices[begin+2])
end


vertices(m::CurvilinearMesh) = m.vertices
faces(m::CurvilinearMesh)    = m.faces

numvertices(m::CurvilinearMesh) = length(m.vertices)
numcells(m::CurvilinearMesh)    = length(m.faces)
# cells(m::CurvilinearMesh)       = Base.OneTo(length(m.faces))
function cells(m::CurvilinearMesh) m.faces end
function cell(m::CurvilinearMesh, i::Int) m.faces[i] end

Base.eltype(::Type{CurvilinearMesh{U,D1,N,T,O}}) where {U,D1,N,T,O} = SVector{N,Int}

meshorder(::CurvilinearMesh{U,D1,N,T,O}) where {U,D1,N,T,O} = O

function skeleton_fast(mesh::CurvilinearMesh, dim::Int)
    skeleton_fast(mesh, Val{dim})
end
