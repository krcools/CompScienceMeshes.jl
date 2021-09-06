
# """
# Turn sheet in a two-sided mesh, i.e. a mesh that contains two copies of the
# same geometric entity, one part of the front of the mesh, the other part of
# the back.
# """
# function twosided(edges::AbstractMesh, faces::AbstractMesh)

#     @assert dimension(edges) == 1
#     @assert vertices(edges) == vertices(faces)

#     back = -deepcopy(edges)
#     bnd = boundary(faces)
#     int = submesh(!in(bnd), back)

#     @show length(int)

#     return union(edges, int)
# end


function twosided(sheet::Mesh)

    Mesh(
        sheet.vertices,
        [sheet.faces; flip.(sheet.faces)])
end

function skeleton_twosided(mesh::AbstractMesh, dim::Int)

    @assert dimension(mesh) == 2
    @assert dim == 1 "only construction of edges supported for now"
    meshdim = dimension(mesh)
    @assert 0 <= dim <= meshdim

    if dim == meshdim
        return mesh
    end

    nc = length(mesh)
    T = SVector{dim+1,Int}
    simplices = zeros(T, nc*binomial(meshdim+1,dim+1))

    n = 1
    for (c,cell) in enumerate(mesh)
        for k in 1:3
            simplex = T(cell[k], cell[mod1(k+1,3)])
            if simplex == sort(simplex)
                simplices[n] = simplex
                n += 1
            end
        end
    end

    resize!(simplices,n-1)

    for simplex in simplices
        pos = findall(==(simplex), simplices)
        length(pos) == 2 || continue
        simplices[pos[2]] = flip(simplices[pos[2]])
    end

    Mesh(vertices(mesh), simplices)

end