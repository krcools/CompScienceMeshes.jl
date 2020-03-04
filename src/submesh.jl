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

    sub2sup = zeros(Int, 0)
    # sup2sub = zeros(Int, length(mesh))

    j = 1
    for i in 1:length(kcells)

        kcell = mesh.faces[i]
        if pred(kcell)
            kcells[j] = kcell
            push!(sub2sup, i)
            # sup2sub[i] = j
            j += 1
        end

    end

    kcells = kcells[1:j-1]

    # sub = Mesh(mesh.vertices, kcells)
    return SubMesh(mesh, sub2sup)
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


mutable struct SubMesh{U,D1,T} <: AbstractMesh{U,D1,T}
    supermesh::AbstractMesh{U,D1,T}
    sub2sup::Vector{Int}
    sup2sub::Vector{Int}
    cells::Vector{SVector{D1,Int}}
end

function SubMesh(supermesh, sub2sup)

    sup2sub = zeros(Int, length(supermesh))
    @inbounds for (i,j) in enumerate(sub2sup)
        sup2sub[j] = i
    end

    Cells = cells(supermesh)[sub2sup]
    SubMesh(supermesh, sub2sup, sup2sub, Cells)
end

vertextype(m::SubMesh) = vertextype(m.supermesh)
celltype(m::SubMesh) = celltype(m.supermesh)
coordtype(m::SubMesh) = coordtype(m.supermesh)
dimension(m::SubMesh) = dimension(m.supermesh)
universedimension(m::SubMesh) = universedimension(m.supermesh)

vertexarray(m::SubMesh) = vertexarray(m.supermesh)
function cellarray(m::SubMesh)
    [cells(m.supermesh)[i][j] for i in m.sub2sup, j in 1:dimension(m)+1]
end

vertices(m::SubMesh) = vertices(m.supermesh)
vertices(m::SubMesh, cell) = vertices(m.supermesh, cell)
numvertices(m::SubMesh) = numvertices(m.supermesh)

cells(m::SubMesh) = m.cells #m.supermesh.faces[m.sub2sup]
numcells(m::SubMesh) = length(m.sub2sup)

issubmesh(sub, sup) = (sub == sup)
issubmesh(sub::SubMesh, sup) = (sub.supermesh == sup)
issubmesh(sub::SubMesh, sup::SubMesh) = (sub == sup)

extend(sm, cell) = sm.sub2sup[cell]
restrict(sm, cell) = sm.sup2sub[cell]

chart(mesh::SubMesh, cell) = simplex(vertices(mesh,cell))
