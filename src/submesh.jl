export submesh

function searchsortedfirst(v, x, lo, hi, by, lt)
    lo = lo-1
    hi = hi+1
    @inbounds while lo < hi-1
        m = (lo+hi)>>>1
        if lt(by(v[m]), x)
        #if lt(o, v[m], x)
            lo = m
        else
            hi = m
        end
    end
    return hi
end

function colsearchsortedfirst(A, col)

    n = size(A,2)
    v =  collect(1:n)
    searchsortedfirst(v, col, 1, n, i->A[:,i], lexless)

end


"""
Returns a mesh on the same vertexbuffer as the input mesh. The submesh
will be a mesh of dimension k containing all the k-cells that are in mesh
and that fulfill the predicate pred.

`pred` is a function with signature `pred(cell) -> Bool` returning true if
the simplex is to be added to the submesh under construction.
"""
function submesh(pred, mesh::Mesh, k)

    kcells = cells(mesh, k)

    j = 1
    for i in 1:length(kcells)

        kcell = kcells[i]
        if pred(kcell)
            kcells[j] = kcell
            j += 1
        end

    end

    kcells = kcells[1:j-1]

    #U = size(mesh.vertices, 1)
    #T = eltype(mesh.vertices)

    Mesh(mesh.vertices, kcells)
end

import CollisionDetection
CD = CollisionDetection


function octree{U,D,T}(mesh::Mesh{U,D,T}, kcells::Array{Int,2})

    nverts = size(kcells,1)
    ncells = size(kcells,2)

    points = zeros(T, ncells, U)
    radii = zeros(T, ncells)

    for i in 1 : ncells
        kcell = kcells[:,i]
        verts = mesh.vertices[:,kcell]

        points[:,i] = sum(verts,2) / nverts, 0.0
        for j in 1 : nverts
            radii[i] = max(radii[i], norm(verts[:,j] ))
        end
    end

    CD.Octree(points, radii)
end

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
Create a predicate that tests wheter a k-cell of bigmesh is
included in smallmesh.
"""
function subset_predicate(bigmesh::Mesh, smallmesh::Mesh, k::Int)

    # build an octree with the k-cells of smallmesh
    kcells = cells(smallmesh, k)
    tree = octree(smallmesh, kcells)

    function pred(kcell)

        # build the patch
        patch(bigmesh, kcell)


    end

    return pred
end
