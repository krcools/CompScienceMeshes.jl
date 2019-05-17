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

    T = eltype(eltype(Γ₁.vertices))
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
