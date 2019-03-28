struct NormalMesh{U,D1,T} <: AbstractMesh{U,D1,T}
    vertices::Vector{SVector{U,T}}
    faces::Vector{SVector{D1,Int}}
    normals::Vector{Int}
end
