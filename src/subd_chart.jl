export subd_chart, chart, getelementVertices
mutable struct subd_chart{T}
    index::Int64
    N::Int64
    RingNodes::Vector{Int64}
    vertices::Vector{SVector{3,T}}
end

domain(ch::subd_chart{T}) where {T}= ReferenceSimplex{2,T,3}()

function chart(Smesh::subdMesh,E)
    T = coordtype(Smesh)
    element = Smesh.elements[E]
    Svertices = Smesh.vertices
    nodes = element.RingNodes
    N = length(nodes)
    verticecoords = Vector{SVector{3,T}}(undef,N)

    for i = 1:N
        verticecoords[i] = Svertices[nodes[i]].Coords
    end
    return chart = subd_chart(E,N,nodes,verticecoords)
end

function getelementVertices(chart::subd_chart{T}) where {T}
    verticescoords = chart.vertices
    verticecoords3 = Vector{SVector{3,T}}(3)
    verticecoords3[1] = verticescoords[1]
    verticecoords3[2] = verticescoords[2]
    verticecoords3[3] = verticescoords[chart.N-6]
    return verticecoords3
end
verticeslist(chart::subd_chart) = chart.vertices