# Mesh: 1D elements embedded in R^U, with N control points per element (N = order+1)
mutable struct CurvilinearMesh{U,N,T,O} <: AbstractMesh{U,N,T}
    vertices::Vector{SVector{U,T}}   # coordinates in R^U
    faces::Vector{SVector{N,Int}}    # connectivity (local order = Gmsh order)

    function CurvilinearMesh(
        vertices::Vector{SVector{U,T}},
        faces::Vector{SVector{N,Int}},
        order::Integer
    ) where {U,N,T}
        order > 0 || throw(ArgumentError("order must be positive."))
        new{U,N,T,Int(order)}(vertices, faces)
    end
end



# Constructors from abstract containers
CurvilinearMesh(vertices::AbstractVector{<:SVector{U,T}},
                faces::AbstractVector{<:SVector{N,Int}},
                order::Integer) where {U,N,T} =
    CurvilinearMesh(Vector{SVector{U,T}}(vertices),
                    Vector{SVector{N,Int}}(faces),
                    order)

# From matrices (verts: U×nV, faces: N×nF)
function CurvilinearMesh(verts::AbstractMatrix{T},
                         faces::AbstractMatrix{<:Integer},
                         order::Integer) where {T}
    U = size(verts, 1)
    N = size(faces, 1)
    v = [SVector{U,T}(verts[:, i]) for i in 1:size(verts, 2)]
    f = [SVector{N,Int}(faces[:, i]) for i in 1:size(faces, 2)]
    CurvilinearMesh(v, f, order)
end

# Mesh interface
coordtype(::CurvilinearMesh{U,N,T,O}) where {U,N,T,O} = T
vertextype(::CurvilinearMesh{U,N,T,O}) where {U,N,T,O} = SVector{U,T}
universedimension(::CurvilinearMesh{U}) where {U} = U
dimension(::CurvilinearMesh) = 1

function celltype(m::CurvilinearMesh{U,N,T,O}) where {U,N,T,O} SimplexGraph{N-O+1} end
function celltype(m::CurvilinearMesh{U,N,T,O}, ::Type{Val{M}}) where {U,N,T,O,M} SimplexGraph{M+1} end

function indextype(m::CurvilinearMesh{U,N}) where {U,N} SVector{N-O+1,Int} end
function indextype(m::CurvilinearMesh{U,N}, ::Type{Val{M}}) where {U,N,M} SVector{M+1,Int} end

function indices(m::CurvilinearMesh{U,N,T,O}, cell) where {U,N,T,O}
    # Currently, we only support lines
    # First come the topological nodes
    return SVector(m.faces[cell][begin], m.faces[cell][begin+1])
end


vertices(m::CurvilinearMesh) = m.vertices
faces(m::CurvilinearMesh)    = m.faces

numvertices(m::CurvilinearMesh) = length(m.vertices)
numcells(m::CurvilinearMesh)    = length(m.faces)
cells(m::CurvilinearMesh)       = Base.OneTo(length(m.faces))
cell(m::CurvilinearMesh, i::Int) = m.faces[i]

Base.eltype(::Type{CurvilinearMesh{U,N,T,O}}) where {U,N,T,O} = SVector{N,Int}

meshorder(::CurvilinearMesh{U,N,T,O}) where {U,N,T,O} = O

function skeleton_fast(mesh::CurvilinearMesh, dim::Int)
    skeleton_fast(mesh, Val{dim})
end
