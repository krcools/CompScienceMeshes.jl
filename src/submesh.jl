export submesh, octree, boundingbox, interior
export overlap_predicate, interior_predicate

using CollisionDetection

"""
Returns a mesh on the same vertexbuffer as the input mesh. The submesh
will be a mesh of dimension k containing all the k-cells that are in mesh
and that fulfill the predicate pred.

`pred` is a function with signature `pred(cell) -> Bool` returning true if
the simplex is to be added to the submesh under construction.
"""
function submesh(pred, mesh::Mesh)

    #kcells = deepcopy(cells(mesh, k))
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
function octree{U,D1,T}(mesh::Mesh{U,D1,T})

    nverts = D1             # number of vertices per cell
    kcells = mesh.faces
    ncells = length(kcells)  # number of cells to store

    points = zeros(Point{U,T}, ncells)
    radii = zeros(T, ncells)

    for i in 1 : ncells
        kcell = kcells[i]
        verts = mesh.vertices[kcell]

        points[i] = sum(verts) / nverts
        for j in 1 : nverts
            radii[i] = max(radii[i], norm(verts[j]-points[i]))
        end
    end

    return Octree(points, radii)
end




"""
Returns the boundingbox of a patch in terms of its center and halfsize.

    function boundingbox{U,D,C,N,T}(p::FlatCellNM{U,D,C,N,T}) -> center, halfsize
"""
function boundingbox{U,D,C,N,T}(p::FlatCellNM{U,D,C,N,T})

    ll = minimum(p.vertices)
    ur = maximum(p.vertices)

    center = (ll + ur) / 2
    halfsize = maximum(ur - center)

    return center, halfsize
end

function boundingbox(V)

    ll = minimum(V)
    ur = maximum(V)

    c = (ll + ur) / 2
    s = maximum(ur - c)

    return c, s
end



"""
Create a predicate that tests wheter a k-cell of bigmesh is
included in smallmesh.
"""
function overlap_predicate(bigmesh::Mesh, smallmesh::Mesh)

    # build an octree with the k-cells of smallmesh
    tree = octree(smallmesh)

    function pred(simplex)

        # create a patch object
        p1 = patch(vertices(bigmesh, simplex), Val{1})

        # find all simplices of the small mesh that potentially
        # collide with this patch
        c1, s1 = boundingbox(p1)
        for box in boxes(tree, (c,s)->boxesoverlap(c,s,c1,s1))
            for i in box.data
                p2 = patch(smallmesh.vertices[smallmesh.faces[i]], Val{1})
                overlap(p1,p2) && return true
            end
        end

        return false
    end

    return pred
end

"""
Creates a predicate that can be used to check wheter an edge is interior to
a surface (true) or on its boundary (false).
"""
function interior_predicate(mesh::Mesh)

    @assert dimension(mesh) == 2

    vtoc, vton = vertextocellmap(mesh)

    function pred(simplex::Vec{2,Int})
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

function interior(mesh::Mesh)
    D = dimension(mesh)
    edges = cells(mesh, D-1)
    pred = interior_predicate(mesh)
    submesh(pred, edges)
end

submesh(sm::Mesh, bm::Mesh) = submesh(overlap_predicate(bm, sm), bm)
