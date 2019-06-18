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

universedimension(m::SComplex2D) = length(vertextype(m))
universedimension(m::SComplex1D) = length(vertextype(m))
universedimension(m::SComplex0D) = length(vertextype(m))

vertextype(::SComplex2D{T,P}) where {T,P} = P
vertextype(::SComplex1D{T,P}) where {T,P} = P
vertextype(::SComplex0D{T,P}) where {T,P} = P

vertices(m::SComplex2D) = m.vertices
vertices(m::SComplex1D) = m.vertices
vertices(m::SComplex0D) = m.vertices
numvertices(m::SComplex2D) = length(m.vertices)

numcells(m::SComplex2D) = length(m.faces)
numcells(m::SComplex1D) = length(m.edges)
numcells(m::SComplex0D) = length(m.nodes)

cells(m::SComplex2D) = m.faces
cells(m::SComplex1D) = m.edges
cells(m::SComplex0D) = m.nodes

celltype(m::SComplex2D) = eltype(m.faces)
celltype(m::SComplex1D) = eltype(m.edges)
celltype(m::SComplex0D) = eltype(m.nodes)

function skeleton(m::SComplex2D, dim::Int)
    T = coordtype(m)
    P = vertextype(m)
    if dim == 2
        return m
    elseif dim == 1
        return SComplex1D{T,P}(m.vertices, m.nodes, m.edges)
    elseif dim == 0
        return SComplex0D{T,P}(m.vertices, m.nodes)
    else
        error("Not a valid value for `dim`")
    end
end

function skeleton(m::SComplex1D, dim::Int)
    T = coordtype(m)
    P = vertextype(m)
    if dim == 1
        return SComplex1D{T,P}(m.vertices, m.nodes, m.edges)
    elseif dim == 0
        return SComplex0D{T,P}(m.vertices, m.nodes)
    else
        error("Not a valid value for `dim`")
    end
end


function chart(m::SComplex2D, face::NTuple{3,Int})

    # Get the comprising vertices in the correct order
    P = vertextype(m)

    # vs = P[]
    vs = Vector{P}(undef,3)
    k = 1
    for i in face
        edge = m.edges[abs(i)]
        if i > 0
            v = m.vertices[m.nodes[edge[1]][1]]
            vs[mod1(k+1,3)] = v
        else
            v = m.vertices[m.nodes[edge[2]][1]]
            vs[mod1(k+1,3)] = v
        end
        k += 1
    end

    @assert k == 4
    return simplex(vs...)
end

function chart(m::SComplex1D, edge::NTuple{2,Int})
    i = edge[1]
    j = edge[2]
    v = m.vertices[m.nodes[i][1]]
    w = m.vertices[m.nodes[j][1]]
    simplex(v,w)
end

function chart(m::SComplex0D, node::NTuple{1,Int})
    i = node[1]
    v = m.vertices[i]
    simplex([v], Val{0})
end


function connectivity(edges::SComplex1D, faces::SComplex2D, op=sign)
    D = spzeros(Int, numcells(faces), numcells(edges))
    for (i,f) in enumerate(cells(faces))
        for (p,j) in enumerate(f)
            D[i,abs(j)] = op(sign(j)*p)
        end
    end
    return D
end


function interior_tpredicate(mesh::SComplex2D)

    edges = skeleton(mesh,1)
    D = connectivity(edges, mesh)
    d = sum(D,dims=1)
    M = Dict((edge,i) for (i,edge) in enumerate(cells(edges)))
    pred(cell) = (d[M[cell]] == 0)
    return pred

end


function boundary(mesh::SComplex2D)

    nodes = skeleton(mesh,0)
    edges = skeleton(mesh,1)

    D = connectivity(edges, mesh)
    d = sum(abs.(D),dims=1)
    M = Dict((edge,i) for (i,edge) in enumerate(cells(edges)))
    pred(cell) = (d[M[cell]] >= 2)
    # pred = interior_tpredicate(mesh)
    Edges = similar(mesh.edges,0)
    for edge in cells(edges)
        !pred(edge) && push!(Edges,edge)
    end

    # Keep only those nodes that are on a retained edge
    kept = falses(numcells(nodes))
    for edge in Edges
        for node in edge
            kept[node] = true
        end
    end

    Nodes = similar(mesh.nodes,0)
    nodes_map = zeros(Int, numcells(nodes))
    num_kept = 0
    for (i,b) in enumerate(kept)
        if b
            num_kept += 1
            nodes_map[i] = num_kept
            push!(Nodes, cells(nodes)[i])
        end
    end

    for (i,edge) in enumerate(Edges)
        Edges[i] = map(p->nodes_map[p], edge)
    end

    T = coordtype(mesh)
    P = vertextype(mesh)
    return SComplex1D{T,P}(mesh.vertices, Nodes, Edges)
end

function CompScienceMeshes.interior(mesh::CompScienceMeshes.SComplex2D, edges=skeleton(mesh,1))
    @assert dimension(mesh) == 2
    @assert vertices(mesh) === vertices(edges)

    C = connectivity(edges, mesh)
    @assert size(C) == (numcells(mesh), numcells(edges))

    nn = vec(sum(abs.(C), dims=1))
    T = CompScienceMeshes.celltype(edges)
    interior_edges = Vector{T}()
    for (i,edge) in pairs(cells(edges))
        nn[i] > 1 && push!(interior_edges, edge)
    end

    # Keep only those nodes that are on a retained edge
    nodes = skeleton(mesh,0)
    kept = falses(numcells(nodes))
    for edge in interior_edges
        for node in edge
            kept[node] = true
        end
    end

    Nodes = similar(mesh.nodes,0)
    nodes_map = zeros(Int, numcells(nodes))
    num_kept = 0
    for (i,b) in enumerate(kept)
        if b
            num_kept += 1
            nodes_map[i] = num_kept
            push!(Nodes, cells(nodes)[i])
        end
    end

    for (i,edge) in enumerate(interior_edges)
        interior_edges[i] = map(p->nodes_map[p], edge)
    end

    T = coordtype(mesh)
    P = vertextype(mesh)
    CompScienceMeshes.SComplex1D{T,P}(mesh.vertices, Nodes, interior_edges)
end

function SComplex2D(mesh::Mesh)
    @assert dimension(mesh) == 2

    mesh0 = skeleton(mesh,0)
    mesh1 = skeleton(mesh,1)

    Nodes = [(c[1],) for c in cells(mesh0)]
    node_map = Dict((n[1],i) for (i,n) in enumerate(Nodes))

    Edges = Tuple{Int,Int}[]
    edge_map = Dict{Tuple{Int,Int},Int}()
    for (i,edge) in enumerate(cells(mesh1))
        n1 = node_map[edge[1]]
        n2 = node_map[edge[2]]
        push!(Edges, (n1,n2,))
        edge_map[(edge[1],edge[2])] = i
        edge_map[(edge[2],edge[1])] = -i
    end

    Faces = Tuple{Int,Int,Int}[]
    for face in cells(mesh)
        v1, v2, v3 = face
        e1 = edge_map[v1,v2]
        e2 = edge_map[v2,v3]
        e3 = edge_map[v3,v1]
        push!(Faces, (e1,e2,e3))
    end

    T, P = coordtype(mesh), vertextype(mesh)
    return SComplex2D{T,P}(mesh.vertices, Nodes, Edges, Faces)
end

function Base.:-(mesh::SComplex2D)
    Faces = similar(mesh.faces,0)
    for f in mesh.faces
        push!(Faces,(-f[1],-f[3],-f[2]),)
    end
    T = coordtype(mesh)
    P = vertextype(mesh)
    return SComplex2D{T,P}(mesh.vertices, mesh.nodes, mesh.edges, Faces)
end
