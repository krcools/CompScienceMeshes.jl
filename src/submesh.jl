export submesh

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
function octree{U,D,T,K1}(mesh::Mesh{U,D,T}, kcells::Array{Vec{K1,Int},1})

    nverts = K1              # number of vertices per cell
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
Predicate used for iteration over an Octree. Returns true if two boxes
specified by their centers and halfsizes overlap. More carefull investigation
of the objects within is required to assess collision.

    function boxesoverlap(c1, hs1, c2, hs2)
"""
function boxesoverlap(c1, hs1, c2, hs2)

    dim = length(c1)
    @assert dim == length(c2)

    hs = hs1 + hs2
    for i in 1 : dim
        if abs(c1[i] - c2[i]) < hs
            return true
        end
    end

    return false
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



"""
Create a predicate that tests wheter a k-cell of bigmesh is
included in smallmesh.
"""
function subset_predicate(bigmesh::Mesh, smallmesh::Mesh)

    # build an octree with the k-cells of smallmesh
    tree = octree(smallmesh, smallmesh.faces)

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


submesh(bm::Mesh, sm::Mesh) = submesh(subset_predicate(bm, sm), bm)
