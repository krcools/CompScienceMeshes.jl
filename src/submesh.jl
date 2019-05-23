using CollisionDetection

"""
Returns a mesh on the same vertexbuffer as the input mesh. The submesh
will be a mesh of dimension k containing all the k-cells that are in mesh
and that fulfill the predicate pred.

`pred` is a function with signature `pred(cell) -> Bool` returning true if
the simplex is to be added to the submesh under construction.
"""
function submesh(pred, mesh::Mesh)

    kcells = similar(mesh.faces)

    j = 1
    for i in 1:length(kcells)

        kcell = mesh.faces[i]
        if pred(kcell)
            kcells[j] = kcell
            j += 1
        end

    end

    kcells = kcells[1:j-1]

    Mesh(mesh.vertices, kcells)
end

"""
Store the k-cells of a mesh in an octree.

    function octree{U,D,T}(mesh::Mesh{U,D,T}, kcells::Array{Int,2})
"""
function octree(mesh)

    nverts = dimension(mesh) + 1 # number of vertices per cell
    ncells = numcells(mesh)

    P = vertextype(mesh)
    T = coordtype(mesh)

    points = zeros(P, ncells)
    radii = zeros(T, ncells)

    for (i,cl) in enumerate(cells(mesh))
        ch = chart(mesh, cl)
        points[i] = cartesian(center(ch))
        radii[i] = maximum(norm(v-points[i]) for v in ch.vertices)
    end

    return Octree(points, radii)
end




"""
Returns the boundingbox of a patch in terms of its center and halfsize.

    function boundingbox{U,D,C,N,T}(p::Simplex{U,D,C,N,T}) -> center, halfsize
"""
function boundingbox(p::Simplex{U,D,C,N,T}) where {U,D,C,N,T}

    # ll = minimum(p.vertices)
    # ur = maximum(p.vertices)

    ll = first(p.vertices); for v in p.vertices; ll = min.(v, ll); end
    ur = first(p.vertices); for v in p.vertices; ll = max.(v, ll); end

    center = (ll + ur) / 2
    halfsize = maximum(ur - center)

    return center, halfsize
end

function boundingbox(V)

    # ll = minimum(V)
    # ur = maximum(V)

    ll = first(V); for v ∈ V; ll = min.(v, ll); end
    ur = first(V); for v ∈ V; ur = max.(v, ur); end

    c = (ll + ur) / 2
    s = maximum(ur - c)

    return c, s
end


boundingbox(p::SVector{N,T}) where {N,T<:Number} = p, zero(eltype(p))


"""
    overlap_gpredicate(γ::Mesh) -> (patch -> Bool)

Create a predicate that for a given patch determinees if it overlaps with the provided target mesh `γ`.
"""
function overlap_gpredicate(γ::Mesh)

    if numcells(γ) == 0
        return simplex -> false
    end

    # build an octree with the k-cells of smallmesh
    tree = octree(γ)

    function pred(p1)

        # find all simplices of the small mesh that potentially
        # collide with this patch
        c1, s1 = boundingbox(p1)
        for box in boxes(tree, (c,s)->boxesoverlap(c,s,c1,s1))
            for i in box
                p2 = chart(γ, γ.faces[i])
                overlap(p1,p2) && return true
            end
        end

        return false
    end

    return pred
end


"""
Geometric predicate for determining in log(N) complexity if a the image of a chart is in
the closure of mesh `γ`.
"""
function inclosure_gpredicate(γ::Mesh)

    if numcells(γ) == 0
        return simplex -> false
    end

    # build an octree with the k-cells of smallmesh
    tree = octree(γ)

    function pred(p1)

        # find all simplices of the small mesh that potentially
        # collide with this patch
        c1, s1 = boundingbox(p1)
        for box in boxes(tree, (c,s)->boxesoverlap(c,s,c1,s1))
            for i in box
                p2 = simplex(vertices(γ, γ.faces[i]))
                isinclosure(p2, p1) && return true
            end
        end

        return false
    end

    return pred
end


"""
Creates a predicate that can be used to check wheter an edge is interior to
a surface (true) or on its boundary (false). This predicate is based on combinatorics. In particular it expects as argument a tuple of indices pointing into the vertex buffer of `mesh`.
"""
function interior_tpredicate(mesh::Mesh)

    @assert dimension(mesh) == 2

    vtoc, vton = vertextocellmap(mesh)

    function pred(simplex::SVector{2,Int})
        v = simplex[1]
        n = vton[v]
        a = vtoc[v,1:n]

        v = simplex[2]
        n = vton[v]
        b = vtoc[v,1:n]

        ab = a ∩ b
        return length(ab) == 2
    end

    return pred
end

"""
Creates a predicate that can be used to check wheter a vertex is interior to
a surface (true) or on its boundary (false).
In particular it expects as argument an index pointing into the vertex buffer of `mesh`.
"""
function interior_vpredicate(mesh::Mesh)
    #TODO: update to accept simplex as argument
    @assert dimension(mesh) == 2
    skel = skeleton(mesh,1)
    edges = cells(skel)
    vtoc,vton = vertextocellmap(skel)
    pred = interior_tpredicate(mesh)
    function vpred(vidx)
        for j = 1:vton[vidx]
            pred(edges[vtoc[vidx,j]]) ? nothing : return false;
        end
        return true
    end
    return vpred
end

"""
    submesh(selection, mesh)

Create a submesh from `mesh` comprising those elements that overlap with elements from `selection`. It is assumed that `selection` and `mesh` have the same dimension.
"""
function submesh(sm::Mesh, bm::Mesh)
    overlaps = overlap_gpredicate(sm)
    in_smallmesh = c -> overlaps(simplex(vertices(bm,c)))
    submesh(in_smallmesh, bm)
end
