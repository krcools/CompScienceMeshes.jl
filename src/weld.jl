using CollisionDetection


"""
    weld(mesh1, mesh2, ...;boundary=false) -> welded_mesh

Build a mesh by welding or pasting together the inputs. Vertices from different
meshes that coincide up to the tolerance will be merged into one. The order
cells appear in the output mesh is equal to the order in the inputs.
boundary = true; will merge only vertices of different meshes if one of
them exists on the 'open' boundary of the mesh.
"""
function weld end


function weld(Γ₁, Γ₂; boundary=false)

    # T = eltype(eltype(Γ₁.vertices))
    T = coordtype(Γ₁)
    tol = sqrt(eps(T))

    verts1 = skeleton(Γ₁,0) # TODO: check this, it might break for unused verts
    radii = zeros(numcells(verts1))
    indcs = [v[1] for v in cells(verts1)]
    cntrs = [cartesian(center(chart(verts1,v))) for v in cells(verts1)]
    tree = Octree(cntrs, radii)

    nv1 = numvertices(Γ₁)
    nv2 = numvertices(Γ₂)

    idmap = collect(nv1 .+ (1:nv2))
    num_equal_vertices = 0

    V1 = vertices(Γ₁)
    V2 = similar(V1,0)

    in_interior1 = interior_vpredicate(Γ₁)
    in_interior2 = interior_vpredicate(Γ₂)
    cpred(u,v,i,j) = boundary ? norm(u-v) < tol && (!in_interior1(indcs[i])
                || !in_interior2(j)) : norm(u-v) < tol

    for (j,v) in enumerate(Γ₂.vertices)
        found = false
        pred(c,s) = fitsinbox(Array(v), 0.0, c, s+tol)
        for box in boxes(tree, pred)
            for i in box
                u = cntrs[i]
                if cpred(u,v,i,j) #norm(u-v) < tol
                    idmap[j] = indcs[i]
                    num_equal_vertices += 1
                    found = true
                    # we implicitly assume that the first copy of a point
                    # appearing in the vertex buffer is the one actually
                    # referred to in the index buffer. DANGEROUS!
                    break
                end
            end
            found && break
        end
        if !found
            idmap[j] -= num_equal_vertices
            push!(V2, v)
        end
    end

    V = [V1; V2]
    F = [Γ₁.faces; Γ₂.faces]

    nc1 = numcells(Γ₁)
    # for (i,c) in enumerate(Γ₂.faces)
    for (i,c) in enumerate(cells(Γ₂))
        F[nc1+i] = map_ids(c, idmap)
    end

    return Mesh(V,unique(F))

end


function weld(G1::SComplex2D, G2::SComplex2D; seam)

    T = coordtype(G1)
    P = vertextype(G1)
    tol = sqrt(eps(T))

    nodes1 = skeleton(G1,0)
    edges1 = skeleton(G1,1)

    node_ctrs1 = [cartesian(center(chart(nodes1,c))) for c in cells(nodes1)]
    edge_ctrs1 = [cartesian(center(chart(edges1,c))) for c in cells(edges1)]

    node_radii1 = zeros(T, length(node_ctrs1))
    edge_radii1 = [volume(chart(edges1,c)) for c in cells(edges1)] / 2

    node_tree1 = Octree(node_ctrs1, node_radii1)
    edge_tree1 = Octree(edge_ctrs1, edge_radii1)

    nodes2 = skeleton(G2,0)
    edges2 = skeleton(G2,1)

    node_is_on_seam = inclosure_gpredicate(seam)
    edge_is_on_seam = overlap_gpredicate(seam)
    # is_interior_edge = interior_tpredicate(G1)
    # is_interior_node = interior_vpredicate(G1)

    Nodes2 = similar(G2.nodes, 0)
    Edges2 = similar(G2.edges, 0)

    num_equal_edges = 0
    edges_map = collect(numcells(edges1) .+ (1:numcells(edges2)))
    for (j,edge2) in enumerate(cells(edges2))
        found = false
        ctr2 = cartesian(center(chart(edges2, edge2)))
        rad2 = volume(chart(edges2,edge2)) / 2
        for box in boxes(edge_tree1, (c,s)->fitsinbox(ctr2, rad2, c, s+tol))
            for i in box
                ctr1 = edge_ctrs1[i]
                # For now we assume edges are the same if their centroids coincide
                if norm(ctr1-ctr2) < tol && edge_is_on_seam(chart(edges1,cells(edges1)[i]))
                # if norm(ctr1-ctr2) < tol && !is_interior_edge(cells(edges1)[i])
                    edges_map[j] = i
                    num_equal_edges += 1
                    found = true
                    break
                end
                found && break
            end
        end
        if !found
            push!(Edges2, edge2)
            edges_map[j] -= num_equal_edges
        else
            @show ctr2[3]
        end
    end

    Faces2 = similar(G2.faces, 0)
    for face in cells(G2)
        push!(Faces2,map(p->sign(p)*getindex(edges_map,abs(p)),face))
    end

    num_equal_nodes = 0
    nodes_map = collect(numcells(nodes1) .+ (1:numcells(nodes2)))
    for (j,node2) in enumerate(cells(nodes2))
        found = false
        ctr2 = cartesian(center(chart(nodes2, node2)))
        rad2 = 0.0
        for box in boxes(node_tree1, (c,s)->fitsinbox(ctr2, rad2, c, s+tol))
            for i in box
                ctr1 = node_ctrs1[i]
                # For now we assume edges are the same if their centroids coincide
                if norm(ctr1-ctr2) < tol && node_is_on_seam(ctr1)
                # if norm(ctr1-ctr2) < tol && !is_interior_node(cells(nodes1)[i][1])
                    nodes_map[j] = i
                    num_equal_nodes += 1
                    found = true
                    break
                end
                found && break
            end
        end
        if !found
            push!(Nodes2, node2)
            nodes_map[j] -= num_equal_nodes
        end
    end

    # Edges2 = similar(G2.edges, 0)
    for (i,edge) in enumerate(Edges2)
        edge_mapped = map(p->sign(p)*getindex(nodes_map,abs(p)),edge)
        Edges2[i] = edge_mapped
        # push!(Edges2, edge_mapped)
    end

    SComplex2D{T,P}(
        vcat(G1.vertices, G2.vertices),
        vcat(G1.nodes, Nodes2),
        vcat(G1.edges, Edges2),
        vcat(G1.faces, Faces2))
end

@generated function map_ids(c::SVector{N,T}, idmap) where {N,T}
    xp = :(SVector{$N,$T}())
    for i in 1:N
        push!(xp.args, :(idmap[c[$i]]))
    end
    return xp
end


function weld(meshes...; boundary=false)
    @assert length(meshes) > 1
    G = weld(meshes[1], meshes[2], boundary=boundary)
    for i in 3:length(meshes)
        G = weld(G, meshes[i], boundary=boundary)
    end
    G
end
