abstract type AbstractRefinement{U,D1,T} <: AbstractMesh{U,D1,T} end

function ==(m1::M, m2::M) where {M<:AbstractRefinement} m1.mesh == m2.mesh && m1.parent == m2.parent end
vertextype(br::AbstractRefinement) = vertextype(br.mesh)
vertices(br::AbstractRefinement) = vertices(br.mesh)
vertices(m::AbstractRefinement, cell) = vertices(m.mesh,cell)
numvertices(br::AbstractRefinement) = numvertices(br.mesh)
cells(br::AbstractRefinement) = cells(br.mesh)
indices(br::AbstractRefinement, p) = indices(br.mesh, p) 
numcells(br::AbstractRefinement) = numcells(br.mesh)
skeleton(br::AbstractRefinement, dim::Int) = skeleton(br.mesh, dim)
chart(br::AbstractRefinement, cell) = chart(br.mesh, cell)
vertexarray(m::AbstractRefinement) = vertexarray(m.mesh)
cellarray(m::AbstractRefinement) = cellarray(m.mesh)
function celltype(m::AbstractRefinement{U,D1}) where {U,D1} celltype(m.mesh) end
function celltype(m::AbstractRefinement{U,D1}, ::Type{Val{M}}) where {U,D1,M} celltype(m.mesh, Val{M}) end
function indextype(m::AbstractRefinement{U,D1}) where {U,D1} indextype(m.mesh) end
function indextype(m::AbstractRefinement{U,D1}, ::Type{Val{M}}) where {U,D1,M} indextype(m.mesh, Val{M}) end


struct RefinementParent{U,D1,T,M,N} <: AbstractMesh{U,D1,T}
    mesh::M
    refinement::N
end

Base.:(==)(m1, m2::RefinementParent) = (m1 == m2.mesh)
Base.:(==)(m1::RefinementParent, m2) = m1.mesh == m2

RefinementParent{U,D1,T}(mesh::M, parent::N) where {U,D1,T,M,N} = RefinementParent{U,D1,T,M,N}(mesh, parent)

parent(m::AbstractRefinement{U,D1,T}) where {U,D1,T} = RefinementParent{U,D1,T}(m.parent, m)
Base.length(m::RefinementParent) = length(m.mesh)
vertices(m::RefinementParent, p) = vertices(m.mesh, p)
cells(m::RefinementParent) = cells(m.mesh)
function chart(m::RefinementParent, ptr) chart(m.mesh, ptr) end

"""
True if m1 is a direct refinement of m2.
"""
function refines(m1::AbstractRefinement, m2::AbstractMesh)
    parent(m1) == nothing && return false
    # return parent(m1) == m2
    m1.parent === m2
end

function refines(m1::AbstractMesh, m2::AbstractMesh)
    return false
end

struct BarycentricRefinement{U,D1,T,M,MP} <: AbstractRefinement{U,D1,T}
    mesh::M
    parent::MP
end

function parent(m::BarycentricRefinement, p)
    numchildren = factorial(dimension(m)+1)
    return mod1(p, numchildren)
end

function children(m::RefinementParent{<:Any,<:Any,<:Any,<:Any,<:BarycentricRefinement}, p)
    numchildren = factorial(dimension(m)+1)
    [(p-1)*numchildren+i for i in 1:numchildren]
end


"""
    barycentric refinement(mesh) -> refined_mesh

Create the mesh obtained by inserting an extra vertex in the barycenters of all cells and
recusively creating fine cells by connecting the barycenter of a k-cell to the already
constructed refined (k-1)-cells on its boundary.
"""
function barycentric_refinement(mesh::AbstractMesh{U,2}; sort=:spacefillingcurve) where U
    T = coordtype(mesh)
    # Get the points and the faces (segments) we will use in the refinement process
    Edges = skeleton(mesh, 1; sort)

    # Note their number and use it to size the new vertex and segment arrays
    NV, NE = numvertices(mesh), numcells(Edges)

    # define the new vertices array that will hold the old points(coarse) and the
    # new points(finer mesh) of vertices
    verts = Array{vertextype(mesh)}(undef, NV+NE)

    # Now we assign the coarse mesh vetrices (old points) to the first half of the
    # array excatly as the original coarse mesh
    for i in 1:NV
        verts[i] = vertices(mesh)[i]
    end

    # Now we create the faces using the new points as well
    # example:
    # oringinal points_index (1,2,3,4) , original faces ((1,2), (2,3), (3,4),(4,1))
    # new points_index(1,2,3,4,5,6,7,8), new_faces((1,5),(5,2),(2,6),(6,3),(3,7),(7,4),(4,8),(8,1))
    #                  a,b      For       E=1 ->   (a,c),(c,b) and so on ..
    edges = Array{celltype(Edges)}(undef, 2NE)
    for (E,p) in enumerate(Edges)
        # the second half of the vertex array holds the refinement points,
        # one at the center of each coarse segment
        verts[NV+E] = cartesian(center(chart(Edges, p)))
        c = NV + E
        a, b = indices(Edges, p)
        edges[2(E-1) + 1] = SimplexGraph(a,c)
        edges[2(E-1) + 2] = SimplexGraph(c,b)
    end
    refmesh = Mesh(verts, edges)

    # return the refinement tagged with its parent so that refines recognises it
    return BarycentricRefinement{U,2,T,typeof(refmesh),typeof(mesh)}(refmesh, mesh)
end

function barycentric_refinement(mesh::AbstractMesh{U,3}; sort=:spacefillingcurve) where U

    D1 = 3
    T = coordtype(mesh)

    edges = skeleton(mesh,1; sort)
    faces = skeleton(mesh,2)

    NV, NE, NF = numvertices(mesh), numcells(edges), numcells(faces)

    nv = NV + NE + NF
    verts = Array{vertextype(mesh)}(undef,nv)
    for V in 1 : numvertices(mesh)
        verts[V] = vertices(mesh)[V]
    end

    # add a vertex in each edge centroid
    # for (E,Edge) in enumerate(cells(edges))
    for (i,E) in enumerate(edges)
        verts[NV+i] = cartesian(center(chart(edges, E)))
    end

    # add a vertex in each face centroid
    # for (F,Face) in enumerate(cells(faces))
    for (i,F) in enumerate(faces)
        verts[NV+NE+F] = cartesian(center(chart(faces, F)))
    end

    D = copy(transpose(connectivity(edges, faces, identity)))
    rows, vals = rowvals(D), nonzeros(D)

    # add six faces in each coarse face
    nf = 6NF
    fcs = zeros(indextype(mesh), nf)
    C = cells(faces)
    for F in 1 : numcells(faces)
        c = NV + NE + F

        for i in nzrange(D, F)
            E = rows[i]
            e = NV + E

            j = abs(vals[i]) # local index of edge E in face F
            a = C[F][mod1(j+1,3)]
            b = C[F][mod1(j+2,3)]

            fcs[6(F-1)+2(j-1)+1] = index(a,e,c)
            fcs[6(F-1)+2(j-1)+2] = index(b,c,e)
        end
    end

    Nodes = skeleton(mesh, 0)
    # node_ctrs = [vertices(Nodes)[node][1] for node in cells(Nodes)]
    node_ctrs = [cartesian(center(chart(Nodes, node))) for node in cells(Nodes)]
    Nodes.faces = Nodes.faces[sort_sfc(node_ctrs)]

    fcs = [SimplexGraph{3}(fc) for fc in fcs]
    fine = Mesh(verts, fcs)
    D = connectivity(Nodes, fine)
    rows, vals = rowvals(D), nonzeros(D)
    sorted_fcs = Vector{celltype(mesh)}()
    for (i,Node) in enumerate(cells(Nodes))
        for k in nzrange(D,i)
            j = rows[k]
            push!(sorted_fcs, fcs[j])
        end
    end

    # sorted_fcs = fcs
    # @show length(fcs)
    # @show length(sorted_fcs)
    @assert length(fcs) == length(sorted_fcs)

    refmesh = Mesh(verts,fcs)
    M = typeof(refmesh)
    N = typeof(mesh)
    # M = AbstractMesh
    # BarycentricRefinement{U,3,T,M}(Mesh(verts, sorted_fcs), mesh)
    BarycentricRefinement{U,3,T,M,N}(refmesh, mesh)
end

function barycentric_refinement(mesh::Mesh{U,4}) where U

    D1 = dimension(mesh)+1
    T = coordtype(mesh)

    edges = skeleton(mesh,1)
    faces = skeleton(mesh,2)
    tetrs = skeleton(mesh,3)

    NV, NE, NF, NT = numvertices(mesh), numcells(edges), numcells(faces), numcells(tetrs)

    nv = NV + NE + NF + NT
    verts = Array{vertextype(mesh)}(undef,nv)
    for V in 1 : numvertices(mesh)
        verts[V] = mesh.vertices[V]
    end

    # add a vertex in each edge centroid
    for (E,Edge) in enumerate(edges)
        verts[NV+E] = cartesian(center(chart(edges, Edge)))
    end

    # add a vertex in each face centroid
    for (F,Face) in enumerate(faces)
        verts[NV+NE+F] = cartesian(center(chart(faces, Face)))
    end

    # add a vertex in each tet centroid
    for (T,Tetr) in enumerate(tetrs)
        verts[NV+NE+NF+T] = cartesian(center(chart(tetrs, Tetr)))
    end

    D = copy(transpose(connectivity(faces, tetrs, identity)))
    C = copy(transpose(connectivity(edges, faces, identity)))
    rows, vals = rowvals(D), nonzeros(D)

    # add 24 tetrs in each coarse tetr
    nt = 24NT
    idx = 1
    fcs = zeros(indextype(mesh), nt)
    for T in 1 : numcells(tetrs)
        c = NV + NE + NF + T # index of the tetr centroid

        # iterate over all faces of tetr T
        for i in nzrange(D, T)
            F = rows[i]
            f = NV + NE + F # f is the vertex index of the face center

            for k in nzrange(C,F)
                E = rowvals(C)[k]
                e = NV + E # e is the vertex index of the edge center

                a = cells(edges)[E][1]
                b = cells(edges)[E][2]

                if sign(nonzeros(C)[k]) < 0
                    a, b = b, a
                end

                fcs[idx+0] = index(a,e,f,c)
                fcs[idx+1] = index(b,f,e,c)

                idx += 2
            end

            # j = abs(vals[i]) # local index of edge E in face F
            # a = mesh.faces[F][mod1(abs(vals[i])+1,3)]
            # b = mesh.faces[F][mod1(abs(vals[i])+2,3)]
            #
            # fcs[6(F-1)+2(j-1)+1] = index(a,e,c)
            # fcs[6(F-1)+2(j-1)+2] = index(b,c,e)
        end
    end

    # Nodes = skeleton(mesh, 0)
    # node_ctrs = [vertices(Nodes)[node][1] for node in cells(Nodes)]
    # Nodes.faces = Nodes.faces[sort_sfc(node_ctrs)]

    fine = Mesh(verts, fcs)
    for (i,fc) in enumerate(fcs)
        # @show fc
        vs = verts[fc]
        t1 = vs[1]-vs[4]
        t2 = vs[2]-vs[4]
        t3 = vs[3]-vs[4]
        if dot(t1 × t2, t3) < 0
            fcs[i] = @SVector[fc[1],fc[2],fc[4],fc[3]]
        end
    end
    
    MR = typeof(fine)
    MP = typeof(mesh)

    BarycentricRefinement{U,4,T,MR,MP}(fine, mesh)
    # D = connectivity(Nodes, fine)
    # rows, vals = rowvals(D), nonzeros(D)
    # sorted_fcs = Vector{celltype(mesh)}()
    # for (i,Node) in enumerate(cells(Nodes))
    #     for k in nzrange(D,i)
    #         j = rows[k]
    #         push!(sorted_fcs, fcs[j])
    #     end
    # end
    #
    # # sorted_fcs = fcs
    # @show length(fcs)
    # @show length(sorted_fcs)
    # @assert length(fcs) == length(sorted_fcs)
    #
    # M = typeof(mesh)
    # BarycentricRefinement{U,3,T,M}(Mesh(verts, sorted_fcs), mesh)
end

children(mesh::Mesh{3,4}, cell) = [24*(cell-1)+1 : 24*cell]
parent(mesh::BarycentricRefinement{3,4}, cell) = div((cell-1),24)+1
parent(mesh::BarycentricRefinement{U,3} where {U}, cell) = div((cell-1),6)+1
parent(mesh::BarycentricRefinement{U,2} where {U}, cell) = div((cell-1),2)+1


struct BisectingRefinement{U,D1,T,M,MP} <: AbstractRefinement{U,D1,T}
    order::Int
    mesh::M
    parent::MP
end

function parent(m::BisectingRefinement, p)
    order = m.order
    numchildren = (div(order * (order+1), 2) + div((order-1) * order,2))
    return mod1(p, numchildren)
end

function children(m::RefinementParent{<:Any,<:Any,<:Any,<:Any,<:BisectingRefinement}, p)
    order = m.refinement.order
    numchildren = (div(order * (order+1), 2) + div((order-1) * order,2))
    [(p-1)*numchildren+i for i in 1:numchildren]
end

function injections(m::RefinementParent{U,D1,T,<:Any,<:BisectingRefinement}, p) where {U,D1,T}

    order = m.refinement.order
    numchildren = (div(order * (order+1), 2) + div((order-1) * order,2))

    h = 1 / order
    # set1 = [
    #     simplex(point(i*h, j*h), point((i+1)*h, j*h), point(i*h, (j+1)*h))
    #     for j in 0:order-1 for i in 0:order-1-j
    # ]
    # set2 = [
    #     simplex(point((i+1)*h, j*h), point((i+1)*h, (j+1)*h), point(i*h, (j+1)*h))
    #      for j in 0:order-2 for i in 0:order-2-j
    # ]
    # set = [set1; set2]

    
    S = Simplex{2,D1-1,2-D1+1,D1,T}
    fcs = Vector{S}()
    for m in 0:order-1
        for n in 0:order-1-m
            push!(fcs, simplex(point(T,m*h,n*h), point(T,(m+1)*h, n*h), point(T, m*h, (n+1)*h)))
            if m + n <= order - 2
                push!(fcs, simplex(point(T,(m+1)*h, n*h), point(T,(m+1)*h, (n+1)*h), point(T, m*h, (n+1)*h)))
            end
        end
    end

    # @show length(fcs)
    # @show numchildren
    # @assert length(fcs) == numchildren
    return fcs
end

"""
    bisecting_refinement(mesh::AbstractMesh{U,3}; order=2) where U


    bisecting_refinement(mesh) -> refinement

Construct a refinement of `mesh` by inserting a new vertex on every existing edge.
Every face is subdived in four small faces.

Only defined for 2D meshes.
"""
function bisecting_refinement(mesh::AbstractMesh{U,3,T}; order=2) where {U,T}

    MP = typeof(mesh)

    edges = skeleton(mesh,1)
    # faces = skeleton(mesh,2)
    faces = mesh

    NV, NE, NF = numvertices(mesh), numcells(edges), numcells(faces)

    verts = deepcopy(vertices(mesh))

    for (i, E) in enumerate(edges)
        ch = chart(edges, E)
        for n in 1:order-1
            p = neighborhood(ch, SVector(1 - n/order))
            push!(verts, cartesian(p))
        end
    end

    idx = Dict{Tuple{Int,Int}, Int}()
    k = 1
    for n in 1:order-2
        for m in 1:order-1-n
            idx[(n,m)] = k
            k += 1
        end
    end

    for F in faces
        ch = chart(faces, F)
        for n in 1:order-2
            for m in 1:order-1-n
                p = neighborhood(ch, SVector(n/order,m/order))
                push!(verts, cartesian(p))
            end
        end
    end

    # Build a matrix that given a coarse face gives
    # access to the coarse edges making up its boundary
    D = connectivity(faces, edges, identity)
    rows, vals = rowvals(D), nonzeros(D)

    # add four faces in each coarse face
    # nf = (div(order * (order+1), 2) + div((order-1) * order,2)) * NF
    S = SimplexGraph{3}
    fcs = Vector{S}()

    for F in 1 : numcells(faces)

        # v = zeros(Int, order+1, order+1)
        v = fill(-200, order+1, order+1)

        # retrieve the indices of the fine vertices in the nodes of face F
        v[order+1,1] = cells(faces)[F][1]
        v[1,order+1] = cells(faces)[F][2]
        v[1,1] = cells(faces)[F][3]

        # retrieve the indices of the fine vertices along the edges of face F
        # e = zeros(Int,3, order-1)
        for i in nzrange(D, F)
            E = rows[i]
            j, σ = abs(vals[i]), sign(vals[i]) # j = +/-1, +/-2, +/-3
            if j == 1
                for n in 1:order-1
                    v[1, order-n+1] = (σ >0) ?
                        NV + (E-1) * (order-1) + n :
                        NV + (E-1) * (order-1) + (order-n)
                end
            elseif j == 2
                for n in 1:order-1
                    v[n+1, 1] = (σ >0) ?
                        NV + (E-1) * (order-1) + n :
                        NV + (E-1) * (order-1) + (order-n)
                end
            elseif j == 3
                for n in 1:order-1
                    v[order-n+1, n+1] = (σ >0) ?
                        NV + (E-1) * (order-1) + n :
                        NV + (E-1) * (order-1) + (order-n)
                end
            end
        end
        # @show e

        # retrieve the indices of the fine vertices in the interior of face F
        for n in 1:order-2
            for m in 1:order-1-n
                k = idx[(n,m)]
                v[n+1,m+1] = NV + (NE * (order-1)) + (F-1) * div((order-1)*(order-2),2) + k
            end
        end

        # @show v

        # Build the fine faces in the coarse face F
        for n in 1:order
            for m in 1:order-n+1
                push!(fcs, S(SVector(v[n,m], v[n+1,m], v[n,m+1])))
                if m + n <= order
                    push!(fcs, S(SVector(v[n+1,m], v[n+1,m+1], v[n,m+1])))
                end
            end
        end
    end

    fine = Mesh(verts, fcs)
    MR = typeof(fine)
    BisectingRefinement{U,3,T,MR,MP}(order, fine, mesh)
end

"""
    los_triangle_center(vertices) -> center

Compute the center for a "line-of-sight" refinement.
"""
function los_triangle_center(r)

    n = (r[1]-r[3]) × (r[2]-r[3])
    n /= norm(n)

    l = [
        norm(r[2]-r[3]),
        norm(r[3]-r[1]),
        norm(r[1]-r[2])]

    α = [
        acos((l[2]^2 + l[3]^2 - l[1]^2)/(2*l[2]*l[3])),
        acos((l[3]^2 + l[1]^2 - l[2]^2)/(2*l[3]*l[1])),
        acos((l[1]^2 + l[2]^2 - l[3]^2)/(2*l[1]*l[2]))]

    if α[1] ≥ π/2
        β = 2/3 * sin(α[2]) * sin(α[3]) / cos(α[2] - α[3])
        t = -β * dot(r[2]-r[3], r[3]-r[1]) / dot(n×(r[2]-r[1]), r[3]-r[1])
        r = r[1] + β*(r[2]-r[1]) + t*cross(n, r[2]-r[1])
    elseif α[2] ≥ π/2
        β = 2/3 * sin(α[3]) * sin(α[1]) / cos(α[3] - α[1])
        t = -β * dot(r[3]-r[1], r[1]-r[2]) / dot(n×(r[3]-r[2]), r[1]-r[2])
        r = r[2] + β*(r[3]-r[2]) + t*cross(n, r[3]-r[2])
    elseif α[3] ≥ π/2
        β = 2/3 * sin(α[1]) * sin(α[2]) / cos(α[1] - α[2])
        t = -β * dot(r[1]-r[2], r[2]-r[3]) / dot(n×(r[1]-r[3]), r[2]-r[3])
        r = r[3] + β*(r[1]-r[3]) + t*cross(n, r[1]-r[3])
    else
        r = (r[1]+r[2]+r[3]) / 3
    end
    return r
end


function lineofsight_refinement(mesh::Mesh{U,3}) where U

    D1 = 3

    edges = skeleton(mesh,1)
    faces = skeleton(mesh,2)

    NV, NE, NF = numvertices(mesh), numcells(edges), numcells(faces)

    nv = NV + NE + NF
    verts = Array{vertextype(mesh)}(undef,nv)
    for V in 1 : numvertices(mesh)
        verts[V] = mesh.vertices[V]
    end

    # add a vertex in each face los-centre
    for (F,Face) in enumerate(cells(faces))
        R = vertices(mesh)[Face]
        # @show F
        losc = los_triangle_center(R)
        verts[NV+NE+F] = losc
    end

    # For each edge, find the neighboring triangles, and connect
    # them with a straight line. At a new vertex at the edge crossing
    C = connectivity(edges, faces, identity)
    rows = rowvals(C)
    # vals = nonzeros(EF)
    for (E,Edge) in  enumerate(cells(edges))
        K = nzrange(C,E)
        # Breaks on open meshes
        @assert length(K) == 2
        Fs = rows[K]
        c1 = verts[NV+NE+Fs[1]]
        c2 = verts[NV+NE+Fs[2]]
        e1 = verts[Edge[1]]
        e2 = verts[Edge[2]]
        c1c2 = c2 - c1
        e1e2 = e2 - e1
        c1e1 = e1 - c1
        l12 = norm(c1c2)
        S = c1c2 * c1c2'
        A = one(S) - S / l12^2
        t = -dot(c1e1, A*e1e2) / dot(e1e2, A*e1e2)
        m = e1 + t * e1e2
        verts[NV+E] = m
    end

    # # add a vertex in each edge centroid
    # for (E,Edge) in enumerate(cells(edges))
    #     verts[NV+E] = cartesian(center(chart(edges, Edge)))
    # end

    D = copy(transpose(C))
    rows = rowvals(D)
    vals = nonzeros(D)
    # vals = nonzeros(D)

    # add six faces in each coarse face
    nf = 6NF
    fcs = zeros(indextype(mesh), nf)
    for F in 1 : numcells(faces)
        c = NV + NE + F

        for i in nzrange(D, F)
            E = rows[i]
            e = NV + E

            j = abs(vals[i]) # local index of edge E in face F
            a = mesh.faces[F][mod1(j+1,3)]
            b = mesh.faces[F][mod1(j+2,3)]

            fcs[6(F-1)+2(j-1)+1] = index(a,e,c)
            fcs[6(F-1)+2(j-1)+2] = index(b,c,e)
        end
    end

    D1 = 3
    T = coordtype(mesh)
    fine = Mesh(verts, fcs)
    M = typeof(fine)
    MP = typeof(mesh)
    BarycentricRefinement{U,D1,T,M,MP}(Mesh(verts, fcs), mesh)
end
