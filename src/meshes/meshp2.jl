struct Mesh2d2p{T,P,G} <: AbstractMesh{3,3,T}
    vertices::Vector{P}
    faces::Vector{G}
end

struct Mesh1d2p{T,P,G} <: AbstractMesh{3,2,T}
    vertices::Vector{P}
    faces::Vector{G}
end

