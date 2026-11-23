"""
    barycentric_refinement(mesh::CurvilinearMesh{U,3}; sort=:spacefillingcurve)

Barycentric refinement of a mesh of curvilinear triangles of order `O`: every triangle is split
into the six triangles spanned by one of its corners, the midpoint of an adjacent edge, and its
centroid, all taken in the reference (barycentric) coordinates of the parent.

The children are curvilinear triangles of the same order `O`, whose interpolation nodes are the
chart of the parent evaluated at the images of the reference nodes. Since the reference triangle
of a child is an affine image of a part of the parent's, the parent's map of degree `O`
restricted to it is again of degree `O`: the children reproduce the curved surface exactly, and
not only at their vertices, as a flat refinement would.

Nodes shared by several children, of the same or of neighbouring parents, are merged, so that the
refined mesh is conforming. They are identified by their exact (rational) barycentric position on
the parent vertex, edge, or cell they lie on, and the parent's own nodes are reused where they
coincide (for `O ≥ 2`, e.g., the midpoints of the edges). The cells are numbered as in the
refinement of flat meshes: the children of the parent cell `p` are `6(p-1)+1:6p`, the pair
`6(p-1)+2j-1`, `6(p-1)+2j` being adjacent to the local edge `j`, opposite to the corner `j`, and
they have the orientation of their parent. The keyword `sort` is accepted for compatibility
with the method for flat meshes and has no effect.
"""
function barycentric_refinement(mesh::CurvilinearMesh{U,3,N,T,O}; sort=:spacefillingcurve) where {U,N,T,O}

    tuples = gmsh_triangle_index_to_tuple[O]
    verts = copy(vertices(mesh))

    # the key of the point with the exact barycentric coordinates λ in the parent cell F, whose
    # corners have the global indices c, on the parent entity it lies on: a corner, an edge
    # (oriented by its global vertex indices, so that neighbouring cells agree), or the cell
    nodekey(F, c, λ) =
        if count(!iszero, λ) == 1
            (:vertex, c[findfirst(!iszero, λ)])
        elseif count(!iszero, λ) == 2
            i, j = findall(!iszero, λ)
            c[i] < c[j] ? (:edge, c[i], c[j], λ[i]) : (:edge, c[j], c[i], λ[j])
        else
            (:cell, F, λ[1], λ[2])
        end

    # the nodes of the parent mesh, at their reference positions (i, j, k)/O
    nodeindex = Dict{Any,Int}()
    for (F, cell) in enumerate(cells(mesh))
        c = cell.indices
        for (r, t) in enumerate(tuples)
            get!(nodeindex, nodekey(F, c, t .// O), c[r])
        end
    end

    corners = (SVector(1 // 1, 0 // 1, 0 // 1), SVector(0 // 1, 1 // 1, 0 // 1), SVector(0 // 1, 0 // 1, 1 // 1))
    centroid = SVector(1 // 3, 1 // 3, 1 // 3)

    fcs = Vector{HOSimplexGraph{3,O,N}}(undef, 6 * numcells(mesh))
    for (F, cell) in enumerate(cells(mesh))
        c = cell.indices
        ch = chart(mesh, F)
        for j in 1:3
            a = corners[mod1(j + 1, 3)]
            b = corners[mod1(j + 2, 3)]
            m = (a + b) / 2                                # the midpoint of the local edge j
            # the two children adjacent to the edge j, with the orientation of the parent (as
            # (a, m, c) and (b, c, m) in the refinement of flat meshes)
            for (s, Λ) in enumerate(((a, m, centroid), (b, centroid, m)))
                idcs = map(tuples) do t
                    λ = (t[1] * Λ[1] + t[2] * Λ[2] + t[3] * Λ[3]) / O
                    get!(nodeindex, nodekey(F, c, λ)) do
                        push!(verts, barytocart(ch, SVector(T(λ[1]), T(λ[2]))))
                        length(verts)
                    end
                end
                fcs[6 * (F - 1) + 2 * (j - 1) + s] = HOSimplexGraph{3,O,N}(SVector{N,Int}(idcs))
            end
        end
    end

    refmesh = CurvilinearMesh(verts, fcs)
    return BarycentricRefinement{U,3,T,typeof(refmesh),typeof(mesh)}(refmesh, mesh)
end

"""
    chartinparent(refinement::BarycentricRefinement{U,3}, cell)

The cell `cell` of the barycentric refinement of a triangle mesh as a flat simplex in the
parametric domain of its parent, `parent(refinement, cell)`: its vertices are the parameters
`(u, v)` of its corners in the parent. Since the children are affine images of parts of the
parent's reference triangle, the parameters of a point of the child are mapped to those in the
parent by `cartesian(chartinparent(refinement, cell), u)`, also for curvilinear cells, whose
corners cannot be located in the parent from their coordinates alone. The children are
numbered as in `barycentric_refinement`.
"""
function chartinparent(refinement::BarycentricRefinement{U,3,T}, cell) where {U,T}
    k = cell - 6 * (parent(refinement, cell) - 1)
    j, s = divrem(k - 1, 2) .+ 1
    corners = (SVector{2,T}(1, 0), SVector{2,T}(0, 1), SVector{2,T}(0, 0))
    a, b = corners[mod1(j + 1, 3)], corners[mod1(j + 2, 3)]
    m = (a + b) / 2
    c = SVector{2,T}(1, 1) / 3
    return s == 1 ? simplex(a, m, c) : simplex(b, c, m)
end
