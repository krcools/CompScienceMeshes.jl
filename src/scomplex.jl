struct SComplex0D{T,P} <: AbstractMesh{3,1,T}
    vertices::Vector{P}
    nodes::Vector{NTuple{1,Int}}
end

struct SComplex1D{T,P} <: AbstractMesh{3,2,T}
    vertices::Vector{P}
    nodes::Vector{NTuple{1,Int}}
    edges::Vector{NTuple{2,Int}}
end

struct SComplex2D{T,P} <: AbstractMesh{3,3,T}
    vertices::Vector{P}
    nodes::Vector{NTuple{1,Int}}
    edges::Vector{NTuple{2,Int}}
    faces::Vector{NTuple{3,Int}}
end

vertextype(::SComplex2D{T,P}) where {T,P} = P

vertices(m::SComplex2D) = m.vertices
numvertices(m::SComplex2D) = length(m.vertices)
numcells(m::SComplex2D) = length(m.faces)

cells(m::SComplex2D) = m.faces
skeleton(m::SComplex2D) = SComplex1D(m.vertices, m.nodes, m.faces)

function chart(m::SComplex2D, face::NTuple{3,Int})

    # Get the comprising vertices in the correct order
    P = vertextype(m)

    vs = P[]
    for i in face
        edge = m.edges[abs(i)]
        if i > 0
            v = m.vertices[edge[1][1]]
        else
            v = m.vertices[edge[2][1]]
        end
        push!(vs,v)
    end

    return simplex(vs...)
end
