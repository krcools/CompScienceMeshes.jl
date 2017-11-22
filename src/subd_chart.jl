export subd_chart, chart, getelementVertices
type subd_chart
    index::Int64
    N::Int64
    RingNodes::Vector{Int64}
    vertices::Vector{SVector{3,Float64}}
end

domain(ch::subd_chart) = ReferenceSimplex{2,Float64,3}()

function chart(Smesh::subdMesh,E)
    element = Smesh.elements[E]
    Svertices = Smesh.vertices
    nodes = element.RingNodes
    N = length(nodes)
    verticecoords = Vector{SVector{3,Float64}}(N)

    for i = 1:N
        verticecoords[i] = Svertices[nodes[i]].Coords
    end
    return chart = subd_chart(E,N,nodes,verticecoords)
end

function getelementVertices(chart::subd_chart)
    verticescoords = chart.vertices
    verticecoords3 = Vector{SVector{3,Float64}}(3)
    verticecoords3[1] = verticescoords[1].Coords
    verticecoords3[2] = verticescoords[2].Coords
    verticecoords3[3] = verticescoords[chart.N-6].Coords
    return verticescoords3
end
