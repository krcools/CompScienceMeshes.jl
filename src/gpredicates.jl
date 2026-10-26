"""
Store the k-cells of a mesh in an octree.

    function octree{U,D,T}(mesh::Mesh{U,D,T}, kcells::Array{Int,2})
"""
function octree(mesh::AbstractMesh)

    nverts = dimension(mesh) + 1 # number of vertices per cell
    ncells = numcells(mesh)

    P = vertextype(mesh)
    T = coordtype(mesh)

    points = zeros(P, ncells)
    radii = zeros(T, ncells)

    # for (i,cl) in enumerate(cells(mesh))
    for (i,p) in enumerate(mesh)
        ch = chart(mesh, p)
        points[i] = cartesian(center(ch))
        radii[i] = maximum(norm(v-points[i]) for v in ch.vertices)
    end

    return Octree(points, radii)
end

"""
    overlap_gpredicate(γ::Mesh) -> (patch -> Bool)

Create a predicate that for a given patch determinees if it overlaps with the provided target mesh `γ`.
"""
function overlap_gpredicate(γ::AbstractMesh)

    if numcells(γ) == 0
        return simplex -> false
    end

    # build an octree with the k-cells of smallmesh
    tree = octree(γ)

    function pred(p1)

        # find all simplices of the small mesh that potentially
        # collide with this patch
        c1, s1 = boundingbox(p1)
        return anysearchtree(i -> overlap(p1, chart(γ, i)), tree, (c1, s1))
    end

    return pred
end


"""
Geometric predicate for determining in log(N) complexity if a the image of a chart is in
the closure of mesh `γ`.
"""
function inclosure_gpredicate(γ::AbstractMesh)

    if numcells(γ) == 0
        return simplex -> false
    end

    # build an octree with the k-cells of smallmesh
    tree = octree(γ)

    function pred(p1)

        # find all simplices of the small mesh that potentially
        # collide with this patch
        c1, s1 = boundingbox(p1)
        return anysearchtree(i -> isinclosure(chart(γ, i), p1), tree, (c1, s1))
    end

    return pred
end
