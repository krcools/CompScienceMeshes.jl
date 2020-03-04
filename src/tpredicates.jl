import Base.in
function in(mesh::CompScienceMeshes.AbstractMesh)
    cells_mesh = sort.(mesh)
    function f(cell)
        sort(cell) in cells_mesh
    end
    return f
end

"""
Creates a predicate that can be used to check wheter an edge is interior to
a surface (true) or on its boundary (false). This predicate is based on combinatorics.
In particular it expects as argument a tuple of indices pointing into the vertex
buffer of `mesh`.
"""
function interior_tpredicate(mesh::Mesh{U,3} where {U})

    # @assert dimension(mesh) == 2

    vtoc, vton = vertextocellmap(mesh)

    function pred(simplex::SVector{2,Int} where {N})
        v = simplex[1]
        n = vton[v]
        a = vtoc[v,1:n]

        v = simplex[2]
        n = vton[v]
        b = vtoc[v,1:n]

        ab = a ∩ b
        return length(ab) > 1
    end

    return pred
end


function interior_tpredicate(mesh::Mesh{U,2} where {U})

    # @assert dimension(mesh) == 2

    vtoc, vton = vertextocellmap(mesh)

    function pred(simplex::SVector{1,Int} where {N})
        v = simplex[1]
        n = vton[v]
        a = vtoc[v,1:n]

        return length(a) > 1
    end

    return pred
end

function interior_tpredicate(mesh::Mesh{U,4} where {U})

    bnd = boundary(mesh)
    bnd_cells = Set(sort(c) for c in cells(bnd))
    function pred(cell::SVector{3,Int})
        !(sort(cell) in bnd_cells)
    end

    return pred
end

"""
Creates a predicate that can be used to check wheter a vertex is interior to
a surface (true) or on its boundary (false).
In particular it expects as argument an index pointing into the vertex buffer of `mesh`.
"""
function interior_vpredicate(mesh::AbstractMesh)
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
