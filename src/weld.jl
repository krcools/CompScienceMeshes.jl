export weld

using CollisionDetection

function weld(Γ₁, Γ₂)

    T = eltype(eltype(Γ₁.vertices))
    tol = sqrt(eps(T))

    radii = zeros(numvertices(Γ₁))
    tree = Octree(Γ₁.vertices, radii)

    nv1 = numvertices(Γ₁)
    nv2 = numvertices(Γ₂)

    idmap = collect(nv1 + (1:nv2))
    for (j,v) in enumerate(Γ₂.vertices)
        pred(c,s) = fitsinbox(Array(v), 0.0, c, s)
        for box in boxes(tree, pred)
            for i in box.data
                u = Γ₁.vertices[i]
                if norm(u-v) < tol
                    idmap[j] = i
                end
            end
        end
    end

    V = [Γ₁.vertices; Γ₂.vertices]
    F = [Γ₁.faces; Γ₂.faces]

    nc1 = numcells(Γ₁)
    for (i,c) in enumerate(Γ₂.faces)
        F[nc1+i] = map_ids(c, idmap)
    end

    return Mesh(V,F)

end

@generated function map_ids{N,T}(c::Vec{N,T}, idmap)
    xp = :(Vec{$N,$T}())
    for i in 1:N
        push!(xp.args, :(idmap[c[$i]]))
    end
    return xp
end
