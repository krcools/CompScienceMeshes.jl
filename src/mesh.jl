export Mesh

using Combinatorics
using Compat.Iterators

export vertexarray, cellarray

"""
    `U::Int`: indicating dimension of the embedding space
    `D1::Int`: one plus the manifold dimension
    `T<:AbstractFloat`: the type of a vertex coordinate
"""
abstract type AbstractMesh{U,D1,T} end


mutable struct Mesh{U,D1,T} <: AbstractMesh{U,D1,T}
    vertices::Vector{SVector{U,T}}
    faces::Vector{SVector{D1,Int}}
    """
    maps a face on its index in enumeration
    """
    dict::Dict{SVector{D1,Int},Int}
end

function Mesh(vertices, faces)
    # dict = Dict((f,i) for (i,f) in enumerate(faces))
    dict = Dict{Int,Int}()
    T = eltype(eltype(vertices))
    U = length(eltype(vertices))
    D1 = length(eltype(faces))
    Mesh{U,D1,T}(vertices, faces, dict)
end

function indices(m::Mesh, cell)
    return m.faces[cell]
end

vertexarray(m::Mesh) = [ v[i] for v in m.vertices, i in 1:universedimension(m) ]
cellarray(m::Mesh) = [ k[i] for k in m.faces, i in 1:dimension(m)+1 ]


"""
    mesh(type, mdim, udim=mdim+1)

Returns an empty mesh with `coordtype` equal to `type`, of dimension `mdim`
and embedded in a universe of dimension `udim`
"""
mesh(T, mdim, udim=mdim+1) = Mesh(Pt{udim,T}[], SVector{mdim+1,Int}[])



"""
    vt = vertextype(mesh)

Returns type of the vertices used to define the cells of the mesh.
"""
vertextype(m::Mesh) = eltype(m.vertices)



"""
    celltype(mesh)

Returns the type of the index tuples stored in the mesh.
"""
celltype(m::Mesh) = eltype(m.faces)



"""
    coordtype(mesh)

Returns `eltype(vertextype(mesh))`
"""
coordtype(m::AbstractMesh{U,D1,T}) where {U,D1,T} = T
# coordtype(m::Mesh) = eltype(vertextype(m))



"""
    dim = dimension(mesh)

Returns the dimension of the mesh. Note that this is
the dimension of the cells, not of the surrounding space.
"""
dimension(m::AbstractMesh{U,D1}) where {U,D1} = D1-1
# dimension(m::Mesh) = length(celltype(m)) - 1



"""
    universedimension(mesh)

Returns the dimension of the surrounding space. Equals
the number of coordinates required to describe a vertex.
"""
universedimension(m::Mesh) = length(vertextype(m))



"""
    vertices(mesh)

Returns an indexable iterable to the vertices of the mesh
"""
vertices(m::Mesh) = m.vertices



#     vertices(mesh, i)
#     vertices(mesh, I)
#
# Select one or multiple vertices from the mesh. In the second
# form, only statically sized arrays are allowed to discourage
# memory allocation. The returned vector in that case will also
# be statically typed and of the same size as `I`.
vertices(m::Mesh, i::Number) = m.vertices[i]
cellvertices(m::AbstractMesh, cell) = vertices(m, indices(m, cell))
@generated function vertices(m::Mesh, I::SVector)
    N = length(I)
    xp = :(())
    for i in 1:N
        push!(xp.args, :(m.vertices[I[$i]]))
    end
    :(SVector($xp))
end



"""
    numvertices(mesh)

Returns the number of vertices in the mesh.

*Note*: this is the number of vertices in the vertex buffer and might include floatin vertices
or vertices not appearing in any cell. In other words the following is not necessarily true:

```julia
    numvertices(mesh) == numcells(skeleton(mesh,0))
```
"""
numvertices(m::Mesh) = length(m.vertices)



"""
    numcells(mesh)

Returns the number of cells in the mesh.
"""
numcells(m::AbstractMesh) = length(m.faces)



"""
    cells(mesh)

Return an iterable collection containing the cells making up the mesh.
"""
cells(mesh::Mesh) = mesh.faces
# function cells(mesh) eachindex(mesh.faces) end


Base.IteratorSize(::AbstractMesh) = Base.HasLength()
Base.length(m::AbstractMesh) = length(cells(m))
# Base.iterate(m::AbstractMesh, state=0) = iterate(cells(m), state)
Base.iterate(m::AbstractMesh, state=0) = iterate(eachindex(cells(m)), state)

# """
#     cellvertices(mesh, i)
#
# Return an indexable collection containing the vertices of cell `i`.
# Shorthand for `vertices(mesh, cells(mesh, i))`.
# """
# cellvertices(mesh,i) = vertices(mesh, cells(mesh, i))



"""
    translate(mesh, v)

Creates a new mesh by translating `mesh` over vector `v`
"""
translate(Γ::AbstractMesh, v) = Mesh([w + v for w in vertices(Γ)], deepcopy(cells(Γ)))



"""
    translate!(mesh, v)

Translates `mesh` over vector `v` inplace.
"""
function translate!(Γ::Mesh, v)
    for i in 1:length(Γ.vertices)
        Γ.vertices[i] += v
    end
    Γ
end


"""
    flip(cell)

Change the orientation of a cell by interchanging the first to indices.
"""
@generated function flip(cell)
    # generate `T(cell[2],cell[1],cell[3],...)`, with `T = typeof(cell)`
    N = length(cell)
    if N <= 1
        return :(cell)
    end
    xp = :($cell(cell[2], cell[1]))
    for i in 3:N
        push!(xp.args, :(cell[$i]))
    end
    xp
end

"""
    flipmesh!(mesh)

Change the orientation of a mesh
"""
function flipmesh!(mesh)
    mesh.faces .= flip.(mesh.faces)
    mesh
end


flipmesh(mesh) = flipmesh!(deepcopy(mesh))


export mirrormesh, mirrormesh!

"""
    mirror(vertex, normal, anchor)

Mirror vertex across a plane defined by its normal and a containing point.
"""
function mirror(vertex, normal, anchor)
    h = dot(vertex - anchor, normal)
    vertex - 2*h * normal
end

function mirrormesh!(mesh, normal, anchor)
    #mesh.vertices .= mirror.(mesh.vertices, normal, anchor)
    for i in eachindex(mesh.vertices)
        mesh.vertices[i] = mirror(mesh.vertices[i], normal, anchor)
    end
    mesh
end

mirrormesh(mesh, args...) = mirrormesh!(deepcopy(mesh), args...)

function rotate!(Γ::Mesh{3}, v)
    α = norm(v)
    u = v/α

    cα = cos(α/2)
    sα = sin(α/2)

    p = point(cα, sα*u[1], sα*u[2], sα*u[3])
    q = point(cα, -sα*u[1], -sα*u[2], -sα*u[3])

    P = @SMatrix [
        +p[1] -p[2] -p[3] -p[4];
        +p[2] +p[1] -p[4] +p[3];
        +p[3] +p[4] +p[1] -p[2];
        +p[4] -p[3] +p[2] +p[1]]

    Q = @SMatrix [
        +q[1] -q[2] -q[3] -q[4];
        +q[2] +q[1] +q[4] -q[3];
        +q[3] -q[4] +q[1] +q[2];
        +q[4] +q[3] -q[2] +q[1]]

    for i in 1:numvertices(Γ)
        r = vertices(Γ,i)
        t = point(0, r[1], r[2], r[3])
        s = P * Q * t
        Γ.vertices[i] = point(s[2], s[3], s[4])
    end
    Γ
end

function rotate(Γ::Mesh{3}, v)
    R = deepcopy(Γ)
    rotate!(R,v)
    R
end


"""
    fliporientation(mesh)

Changes the mesh orientation inplace. If non-orientatble, undefined.
"""
function fliporientation!(m::Mesh)
    for i in 1:numcells(m)
        m.faces[i] = fliporientation(m.faces[i])
    end
    return m
end


"""
    fliporientation(mesh)

Returns a mesh of opposite orientation.
"""
function fliporientation(m::Mesh)
    n = deepcopy(m)
    fliporientation!(n)
end

@generated function fliporientation(I::SVector{N,T}) where {N,T}
    @assert N >= 2
    xp = :(SVector{N,T}(I[2],I[1]))
    for i in 3:N
        push!(xp.args, :(I[$i]))
    end
    return xp
end

"""
    -mesh -> flipped_mesh

Create a mesh with opposite orientation.
"""
Base.:-(m::Mesh) = fliporientation(m)


Base.getindex(m::AbstractMesh, I::Vector{Int}) = Mesh(vertices(m), cells(m)[I])


"""
    boundary(mesh)

Returns the boundary of `mesh` as a mesh of lower dimension.
"""
function boundary(mesh)

    D = dimension(mesh)

    # vertices have no boundary
    @assert 0 < D

    I = eltype(cells(mesh))
    length(mesh) == 0 && return Mesh(vertices(mesh), I[])

    # build a list of D-1 cells
    edges = skeleton_fast(mesh, D-1)
    faces = skeleton_fast(mesh, D)

    # get the edge-face connection matrix
    conn = connectivity(edges, faces, identity)

    # find the edges that only have one adjacent face
    #i = find(x -> x < 2, sum(abs.(conn), dims=1))
    rows = rowvals(conn)
    vals = nonzeros(conn)

    I = celltype(edges)
    bnd_edges = Vector{I}(undef, length(edges))
    i = 1
    for (e,edge) in enumerate(edges)
        nzr = nzrange(conn,e)
        length(nzr) != 1 && continue
        relop = vals[nzr[1]]
        inds = indices(edges, edge)
        bnd_edges[i] = (relop > 0) ? inds : flip(inds)
        i += 1
    end

    resize!(bnd_edges, i-1)
    bnd = Mesh(vertices(mesh), bnd_edges)
end


"""
Complement to boundary. This function selects those edges that have at
least two faces adjacent. The case with more than two neighboring faces
occurs on non-manifold structures (e.g. containing junctions)
"""
function interior(mesh::Mesh, edges=skeleton(mesh,1))
    @assert dimension(mesh) == 2
    @assert vertices(mesh) === vertices(edges)

    C = connectivity(edges, mesh)
    @assert size(C) == (numcells(mesh), numcells(edges))

    nn = vec(sum(abs.(C), dims=1))
    T = CompScienceMeshes.celltype(edges)
    interior_edges = Vector{T}()
    for (i,edge) in pairs(cells(edges))
        nn[i] > 1 && push!(interior_edges, edge)
    end
    Mesh(vertices(mesh), interior_edges)
end



"""
    vertextocellmap(mesh) -> vertextocells, numneighbors

Computed an V×M array `vertextocells` where V is the number of vertices
and M is the maximum number of cells adjacent to any given vertex such
that `vertextocells[v,i]` is the index in the cells of `mesh` of the `i`th
cell adjacent to teh `v`-th vertex. `numneighbors[v]` contains the number
of cells adjacent to the `v`-th vertex.

This method allows e.g. for the efficient computation of the connectivity
matrix of the mesh.
"""
function vertextocellmap(mesh)

    numverts = numvertices(mesh)
    numcells = length(mesh)
    numneighbors = zeros(Int, numverts)
    for i in mesh
        cell = indices(mesh, i)
        for v in cell
            numneighbors[v] += 1
        end
    end

    npos = -1
    vertstocells = fill(npos, numverts, maximum(numneighbors))
    numneighbors = zeros(Int, numverts)

    for i in mesh
        cell = indices(mesh, i)
        for v in cell
            k = (numneighbors[v] += 1)
            vertstocells[v,k] = i
        end
    end

    vertstocells, numneighbors
end


function vertextocell(mesh)

    row = Dict{Int,Tuple{Int,Int}}()
    k = 1
    for cell in mesh
        inds = indices(mesh, cell)
        for v in inds
            if haskey(row,v)
                (i,n) = row[v]
                row[v] = (i,n+1)
            else
                row[v] = (k,1)
                k += 1
            end
        end
    end

    NC = fill(-1, length(row))
    for (i,n) in values(row)
        NC[i] = n
    end
    @assert !any(NC .== - 1)
    VC = fill(-1, length(row), maximum(NC))

    fill!(NC, 0)
    for (c,cell) in enumerate(mesh)
        inds = indices(mesh, cell)
        for v in inds
            i, _ = row[v]
            k = NC[i] + 1
            VC[i,k] = c
            NC[i] = k
        end
    end

    LG = fill(-1, length(row))
    for (v,(i,n)) in row
        LG[i] = v
    end
    @assert !any(LG .== -1)

    GL = Dict{Int,Int}()
    for (v,(i,n)) in row
        GL[v] = i
    end

    return VC, NC, LG, GL
end



# function faces(s::SVector{4,Int})
#     return [
#         @SVector[4,3,2],
#         @SVector[1,3,4],
#         @SVector[1,4,2],
#         @SVector[1,2,3]
#     ]
# end
#
# function faces(s::SVector{3,Int})
#     return [
#         @SVector[2,3],
#         @SVector[3,1],
#         @SVector[1,2],
#     ]
# end
#
# function faces(s::SVector{2,Int})
#     return[
#         @SVector[2],
#         @SVector[1]
#     ]
# end

"""
    skeleton(mesh, dim)

Returns a mesh comprising the `dim`-dimensional sub cells of `mesh`. For example to retrieve
the edges of a given surface `mesh`,

```julia
edges = skelton(mesh, 1)
```
"""
function skeleton(mesh, dim::Int; sort=:spacefillingcurve)

    meshdim = dimension(mesh)
    @assert 0 <= dim <= meshdim

    if dim == meshdim
        return mesh
    end

    sk = skeleton_fast(mesh,dim)
    sort != :spacefillingcurve && return sk
    
    # sort the simplices on a SFC
    simplices = cells(sk)
    ctrs = [sum(cellvertices(sk,c))/(dim+1) for c in sk]
    if length(sk) > 0
        simplices = simplices[sort_sfc(ctrs)]
    end

    Mesh(vertices(mesh), simplices)
end


function skeleton_fast(mesh, dim::Int)

    meshdim = dimension(mesh)
    @assert 0 <= dim <= meshdim

    if dim == meshdim
        return mesh
    end

    nc = numcells(mesh)
    C = SVector{dim+1,Int}
    simplices = zeros(C, nc*binomial(meshdim+1,dim+1))

    n = 1
    for c = 1 : nc

        # cell = cells(mesh)[c]
        cell = indices(mesh, c)
        for simplex in combinations(cell,dim+1)
            simplices[n] = sort(simplex)
            n += 1
        end
    end

    simplices = unique(simplices)
    Mesh(vertices(mesh), simplices)
end




"""
    skeleton(pred, mesh, dim)

Like `skeleton(mesh, dim)`, but only cells for which `pred(cell)`
returns true are withheld.
"""
# function skeleton(pred, mesh, dim)

#     meshdim = dimension(mesh)
#     @assert 0 <= dim <= meshdim

#     nc = numcells(mesh)
#     C = SVector{dim+1,Int}
#     simplices = zeros(C, nc*binomial(meshdim+1,dim+1))

#     n = 1
#     for c = 1 : nc

#         cell = mesh.faces[c]
#         for simplex in combinations(cell,dim+1)
#             if pred(mesh, SVector{dim+1,Int}(simplex))
#                 simplices[n] = sort(simplex)
#                 n += 1
#             end
#         end
#     end

#     simplices = simplices[1:n-1]
#     simplices = unique(simplices)

#     Mesh(mesh.vertices, simplices)
# end







"""
    connectivity(faces, cells, op=sign)

Create a sparse matrix `D` of size `numcells(cells)` by `numcells(faces)` that
contiains the connectivity info of the mesh. In particular `D[m,k]` is `op(r)`
where `r` is the local index of face `k` in cell `m`. The sign of `r` is
positive or negative depending on the relative orientation of face `k` in cell
`m`.

For `op=sign`, the matrix returned is the classic connectivity matrix, i.e.
the graph version of the exterior derivative.
"""
function connectivity(kcells::AbstractMesh, mcells::AbstractMesh, op = sign)

    vtok, _, lgk, glk = vertextocell(kcells)
    vtom, _, lgm, glm = vertextocell(mcells)

    npos = -1

    dimk = numcells(kcells)
    dimm = numcells(mcells)

    Rows = Int[]
    Cols = Int[]
    Vals = Int[]

    sh = 2 * max(dimm, dimk)
    # sizehint!(Rows, sh)
    # sizehint!(Cols, sh)
    # sizehint!(Vals, sh)

    for vk in axes(vtok,1)
        V = lgk[vk]
        haskey(glm,V) || continue
        vm = glm[V]
        for q in axes(vtok,2)
            i = vtok[vk,q]
            i == npos && break
            # kcell = cells(kcells)[i]
            kcell = indices(kcells,i)
            for s in axes(vtom,2)
                j = vtom[vm,s]
                j == npos && break
                # mcell = cells(mcells)[j]
                mcell = indices(mcells, j)
                val = op(relorientation(kcell, mcell))
                iszero(val) && continue
                push!(Rows, j)
                push!(Cols, i)
                push!(Vals, op(relorientation(kcell, mcell)))
            end
        end
    end

    D = sparse(Rows, Cols, Vals, dimm, dimk, (x,y)->y)
    return D
end

# function connectivity2(kcells::AbstractMesh, mcells::AbstractMesh, op = sign)

#     # if dimension(kcells) > dimension(mcells)
#     #     C = connectivity_impl(mcells, kcells, op)
#     #     return copy(transpose(C))
#     # end
#     return connectivity_impl(kcells, mcells, op)
# end

# import InteractiveUtils

function connectivity2(kcells::AbstractMesh, mcells::AbstractMesh)
    
    @assert dimension(kcells) < dimension(mcells)

    npos = -1
    MCells = cells(mcells)
    vtom, _ = vertextocellmap(mcells)

    Rows = Int[]
    Cols = Int[]
    Vals = Int[]

    for (j,kcell) in enumerate(kcells)
        for v in kcell
            for i in vtom[v,:]
                i == npos && break
                mcell = MCells[i]
                issubset(kcell, mcell) || continue

                push!(Rows, i)
                push!(Cols, j)
                σ = relorientation(kcell, mcell)
                push!(Vals, σ)
            end
        end
    end

    sparse(Rows, Cols, Vals, length(mcells), length(kcells))
end




"""
    pairs = cellpairs(mesh, edges, dropjunctionpair=false)

Given a mesh and set of oriented edges from that mesh (as generated by `skeleton`),
    `cellpairs` will generate a 2 x K matrix, where K is the number of pairs
    and each column contains a pair of indices in the cell array of `mesh` that have
    one of the supplied edges in common.

Returns an array of pairs of indices, each pair corresponding to a pair of adjacent faces.

(If the mesh is oriented, the first row of `facepairs` will contain indices to the cell
    for which the corresponding edge has a positive relative orientation.

If a edge lies on the boundary of the mesh, and only has one neighboring cell, the
    second row of `facepairs` will contain `-k` with `k` the local index of the corresponding
    edge in its neighboring triangle.

If an edge has more than two neighboring cells (i.e. the edge is on a junction),
    all possible pairs of cells that have the junction edge in common are supplied. if
    `dropjunctionpair == false` then one of the possible pairs of cells is not recorded.
    This is done to avoid the creation of linearly dependent basis functions in the
    construction of boundary element methods for Maxwell's equations.)
"""
function cellpairs(mesh, edges; dropjunctionpair=false)

    ndrops = dropjunctionpair ? 1 : 0
    @assert dimension(edges)+1 == dimension(mesh)

    numedges = numcells(edges)

    # perform a dry run to determine the number of cellpairs
    v2e, nn = vertextocellmap(mesh)
    k = 0
    for e in edges
        edge = indices(edges, e)

        v = edge[1]; n = nn[v]; nbd1 = v2e[v,1:n]
        v = edge[2]; n = nn[v]; nbd2 = v2e[v,1:n]

        nbd = intersect(nbd1, nbd2)
        n = length(nbd)
        if n < 3
            k += 1
        else
            k += (n-ndrops)
        end

    end
    facepairs = zeros(Int, 2, k)

    k = 1
    Cells = cells(mesh)
    for e in edges
        edge = indices(edges, e)

        # neighborhood of startvertex
        v = edge[1]
        n = nn[v]
        nbd1 = v2e[v,1:n]

        # neighborhood of endvertex
        v = edge[2]
        n = nn[v]
        nbd2 = v2e[v,1:n]

        # get the intersection
        nbd = intersect(nbd1, nbd2)
        n = length(nbd)
        @assert 0 < n

        if n == 1 # boundary edge

            c = nbd[1]
            cell = Cells[c]
            s = relorientation(edge, cell)
            facepairs[1,k] = c
            facepairs[2,k] = -abs(s)

            k += 1
        elseif n == 2

            c = nbd

            cell1 = Cells[c[1]] #mesh.faces[c[1]]
            cell2 = Cells[c[2]] #mesh.faces[c[2]]

            r1 = relorientation(edge, cell1)
            r2 = relorientation(edge, cell2)

            if r1 > 0 && r2 < 0
                c1 = c[2]
                c2 = c[1]
            else
                c1 = c[1]
                c2 = c[2]
            end

            facepairs[1,k] = c1
            facepairs[2,k] = c2

            k += 1
        else

            for c in drop(combinations(nbd,2), ndrops)

                cell1 = Cells[c[1]] #mesh.faces[c[1]]
                cell2 = Cells[c[2]] #mesh.faces[c[2]]

                r1 = relorientation(edge, cell1)
                r2 = relorientation(edge, cell2)

                if r1 > 0 && r2 < 0
                    c1 = c[2]
                    c2 = c[1]
                else
                    c1 = c[1]
                    c2 = c[2]
                end

                facepairs[1,k] = c1
                facepairs[2,k] = c2

                k += 1
            end
        end
    end

    facepairs
end


"""
    chart(mesh, cell) -> cell_chart

Return a chart describing the supplied cell of `mesh`.
"""
chart(mesh::Mesh, cell) = simplex(vertices(mesh, indices(mesh,cell)))

parent(mesh::AbstractMesh) = nothing
# parent(mesh::Mesh) = nothing

"""
    isoriented(mesh) -> Bool

Returns true is all cells are consistently oriented, false otherwise.
"""
function isoriented(m::AbstractMesh)

    @assert dimension(m) >= 0
    edges = skeleton(m, dimension(m)-1)

    D = connectivity(edges, m)
    S = (abs.(sum(D,dims=1)) .<= 1)
    return all(S)
end

"""
True if m1 is a direct refinement of m2.
"""
function refines(m1::AbstractMesh, m2::AbstractMesh)
    parent(m1) == nothing && return false
    return parent(m1) == m2
end


"""
    union(mesh1, mesh2, ...)

Create the topological union of two meshes. This requires them to be
defined on the same vertex set. No geometric considerations are taken
into account.
"""
function Base.union(m1::AbstractMesh, m2::AbstractMesh)
    @assert dimension(m1) == dimension(m2)
    @assert vertices(m1) == vertices(m2)

    Verts = vertices(m1)
    Cells1 = cells(m1)
    Cells2 = cells(m2)
    Mesh(Verts, vcat(Cells1, Cells2))
end

function Base.union(m1::AbstractMesh, ms::Vararg{AbstractMesh})
    return union(m1, union(ms...))
end



function breadthfirst(op, mesh, root)

    dim = dimension(mesh)
    edges = skeleton(mesh, dim-1)
    D = connectivity(edges, mesh)
    A = D*D'

    rows = SparseArrays.rowvals(A)
    vals = SparseArrays.nonzeros(A)

    discovered = falses(length(mesh))
    queue = DataStructures.Queue{typeof(root)}()
    facing = DataStructures.Queue{Symbol}()

    DataStructures.enqueue!(queue,root)
    DataStructures.enqueue!(facing,:up)
    discovered[root] = true

    while !isempty(queue)

        x = DataStructures.dequeue!(queue)
        s = DataStructures.dequeue!(facing)
        op(mesh, x, s)
        for k in SparseArrays.nzrange(A,x)
            r = rows[k]
            v = vals[k]
            r == x && continue
            discovered[r] && continue
            if v != -1
                q = (s == :up) ? :down : :up
            else
                q = (s == :up) ? :up : :down
            end
            DataStructures.enqueue!(queue,r)
            DataStructures.enqueue!(facing,q)
            discovered[r] = true
        end
    end

    @assert all(discovered)
end


function orient(mesh::Mesh)

    root = 1
    C = collect(cells(mesh))
    function f(mesh,i,v)
        v == :up && return
        cell = C[i]
        mesh.faces[i] = fliporientation(cell) 
    end

    breadthfirst(f, mesh, root)

end

import Base.convert

"""
    convert(::Type{NewT}, mesh::Mesh{U,D1,T}) where {NewT<:Real,U,D1,T}

Converts a Mesh with coordtype T to a Mesh with coordtype NewT.
"""
function convert(::Type{NewT}, mesh::Mesh{U,D1,T}) where {NewT<:Real,U,D1,T}
    vertices=[convert(SVector{U,NewT}, vertex) for vertex in mesh.vertices]
    return Mesh(vertices, mesh.faces)
end
