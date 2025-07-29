using CollisionDetection

"""
Returns a mesh on the same vertexbuffer as the input mesh. The submesh
will be a mesh of dimension k containing all the k-cells that are in mesh
and that fulfill the predicate pred.

`pred` is a function with signature `pred(cell) -> Bool` returning true if
the simplex is to be added to the submesh under construction.
"""
function submesh(pred, mesh::AbstractMesh)

    kcells = cells(mesh)
    sub2sup = zeros(Int, 0)

    j = 0
    # for (i,kcell) = enumerate(cells(mesh))
    for i in mesh
        # kcell = indices(mesh, i)
        if pred(mesh, i)
            push!(sub2sup, i)
            j += 1
        end

    end

    kcells = kcells[1:j]
    return SubMesh(mesh, sub2sup)
end

"""
    submesh(selection, mesh)

Create a submesh from `mesh` comprising those elements that overlap with elements from `selection`. It is assumed that `selection` and `mesh` have the same dimension.
"""
function submesh(sm::Mesh, bm::Mesh)
    overlaps = overlap_gpredicate(sm)
    # in_smallmesh = c -> overlaps(simplex(vertices(bm,c)))
    in_smallmesh = (m,c) -> overlaps(chart(m,c))
    submesh(in_smallmesh, bm)
end


mutable struct SubMesh{U,D1,T} <: AbstractMesh{U,D1,T}
    supermesh::AbstractMesh{U,D1,T}
    sub2sup::Vector{Int}
    sup2sub::Vector{Int}
    cells::Vector{SVector{D1,Int}}
end

# function celltype(m::SubMesh{U,D1}) where {U,D1} SimplexGraph{D1} end
function celltype(m::SubMesh{U,D1}, ::Type{Val{M}}) where {U,D1,M} SimplexGraph{M+1} end
function indextype(m::SubMesh{U,D1}) where {U,D1} SVector{D1,Int} end
function indextype(m::SubMesh{U,D1}, ::Type{Val{M}}) where {U,D1,M} SVector{M+1,Int} end


function SubMesh(supermesh, sub2sup)

    sup2sub = zeros(Int, length(supermesh))
    @inbounds for (i,j) in enumerate(sub2sup)
        sup2sub[j] = i
    end

    # Cells = cells(supermesh)[sub2sup]
    Cells = [indices(supermesh, i) for i in supermesh]
    Cells = Cells[sub2sup]
    SubMesh(supermesh, sub2sup, sup2sub, Cells)
end

function vertextype(m::SubMesh) vertextype(m.supermesh) end
function celltype(m::SubMesh) celltype(m.supermesh) end
function coordtype(m::SubMesh) coordtype(m.supermesh) end 
function dimension(m::SubMesh) dimension(m.supermesh) end
function universedimension(m::SubMesh) universedimension(m.supermesh) end
function parent(m::SubMesh) m.supermesh end

vertexarray(m::SubMesh) = vertexarray(m.supermesh)
function cellarray(m::SubMesh)
    [cells(m.supermesh)[i][j] for i in m.sub2sup, j in 1:dimension(m)+1]
end

vertices(m::SubMesh) = vertices(m.supermesh)
vertices(m::SubMesh, cell) = vertices(m.supermesh, cell)
numvertices(m::SubMesh) = numvertices(m.supermesh)

indices(m::SubMesh, p) = indices(m.supermesh, m.sub2sup[p])
cells(m::SubMesh) = m.cells #m.supermesh.faces[m.sub2sup]
numcells(m::SubMesh) = length(m.sub2sup)

issubmesh(sub, sup) = (sub == sup)
issubmesh(sub::SubMesh, sup) = (sub.supermesh == sup)
issubmesh(sub::SubMesh, sup::SubMesh) = (sub == sup)

extend(sm, cell) = sm.sub2sup[cell]
restrict(sm, cell) = sm.sup2sub[cell]

# chart(mesh::SubMesh, cell) = simplex(vertices(mesh,cell))
chart(mesh::SubMesh, p) = chart(mesh.supermesh, mesh.sub2sup[p])
