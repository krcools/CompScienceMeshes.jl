mutable struct BarycentricRefinement{U,D1,T,M,MP} <: AbstractMesh{U,D1,T}
    mesh::M
    parent::MP
end

vertextype(br::BarycentricRefinement) = vertextype(br.mesh)
vertices(br::BarycentricRefinement) = vertices(br.mesh)
vertices(m::BarycentricRefinement, cell) = vertices(m.mesh,cell)
numvertices(br::BarycentricRefinement) = numvertices(br.mesh)
cells(br::BarycentricRefinement) = cells(br.mesh)
indices(br::BarycentricRefinement, p) = indices(br.mesh, p) 
numcells(br::BarycentricRefinement) = numcells(br.mesh)
skeleton(br::BarycentricRefinement, dim::Int) = skeleton(br.mesh, dim)
chart(br::BarycentricRefinement, cell) = chart(br.mesh, cell)
vertexarray(m::BarycentricRefinement) = vertexarray(m.mesh)
cellarray(m::BarycentricRefinement) = cellarray(m.mesh)

function parent(m::BarycentricRefinement, p)
    numchildren = factorial(dimension(m)+1)
    return mod1(p, numchildren)
end

struct BarycentricRefinementParent{U,D1,T,M,N} <: AbstractMesh{U,D1,T}
    mesh::M
    refinement::N
end

Base.:(==)(m1, m2::BarycentricRefinementParent) = (m1 == m2.mesh)
Base.:(==)(m1::BarycentricRefinementParent, m2) = m1.mesh == m2

BarycentricRefinementParent{U,D1,T}(mesh::M, parent::N) where {U,D1,T,M,N} = BarycentricRefinementParent{U,D1,T,M,N}(mesh, parent)

parent(m::BarycentricRefinement{U,D1,T}) where {U,D1,T} = BarycentricRefinementParent{U,D1,T}(m.parent, m.mesh)
Base.length(m::BarycentricRefinementParent) = length(m.mesh)
vertices(m::BarycentricRefinementParent, p) = vertices(m.mesh, p)
cells(m::BarycentricRefinementParent) = cells(m.mesh)


function children(m::BarycentricRefinementParent, p)
    numchildren = factorial(dimension(m)+1)
    [(p-1)*numchildren+i for i in 1:numchildren]
end


"""
    barycentric refinement(mesh) -> refined_mesh

Create the mesh obtained by inserting an extra vertex in the barycenters of all cells and
recusively creating fine cells by connecting the barycenter of a k-cell to the already
constructed refined (k-1)-cells on its boundary.
"""
function barycentric_refinement(mesh::Mesh{U,2}; sort=:spacefillingcurve) where U

    # Get the points and the faces (segments) we will use in the refinement process
    Edges = skeleton(mesh, 1; sort)

    # Note thier number and use it to find the new number of points and segments(faces)
    # after refienmnts
    NV, NE = numvertices(Edges), numcells(Edges)

    nv = NV + NE
    ne = 2NE

    # define the new vertices array that will hold the old points(coarse) and the
    # new points(finer mesh) of vertices
    verts = similar(mesh.vertices, nv)

    # Now we assign the coarse mesh vetrices (old points) in the upper half excatly
    # as the original coarse mesh
    for i in 1:NV
        verts[i] = mesh.vertices[i]
    end

    # for the second half we store the refienment points and we calculate them using
    # the old points. We insert a point between every pair of points in the original
    # coarse mesh by (point1+ point2)\2
    # for (i, Face) in enumerate(cells(mesh))
    for (i,f) in enumerate(mesh)
        verts[NV+i] = cartesian(center(chart(mesh, f)))
    end

    # Now we create the faces using the new points as well
    # example:
    # oringinal points_index (1,2,3,4) , original faces ((1,2), (2,3), (3,4),(4,5))
    # new points_index(1,2,3,4,5,6,7,8), new_faces((1,5),(5,2),(2,6),(6,3),(3,7),(7,4),(4,8),(8,1))
    #                  a,b      For       E=1 ->   (a,c),(c,b) and so on ..
    edges = similar(Edges.faces, ne)
    for E in 1:NE
        c = NV + E
        a = Edges.faces[E][1]
        b = Edges.faces[E][2]
        edges[2(E-1) + 1] = SVector(a,c)
        edges[2(E-1) + 2] = SVector(c,b)
    end

    # return the new mesh after refienments
    return Mesh(verts, edges)
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
    fcs = zeros(celltype(mesh), nf)
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
    node_ctrs = [vertices(Nodes)[node][1] for node in cells(Nodes)]
    Nodes.faces = Nodes.faces[sort_sfc(node_ctrs)]

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
    fcs = zeros(celltype(mesh), nt)
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

"""
    bisecting_refinement(mesh) -> refinement

Construct a refinement of `mesh` by inserting a new vertex on every existing edge.
Every face is subdived in four small faces.

Only defined for 2D meshes.
"""
function bisecting_refinement(mesh::Mesh{U,3}) where U

    D1 = 3

    edges = skeleton(mesh,1)
    faces = skeleton(mesh,2)

    NV, NE, NF = numvertices(mesh), numcells(edges), numcells(faces)

    nv = NV + NE
    verts = Array{vertextype(mesh)}(undef,nv)

    # copy over the existing vertices
    for V in 1 : numvertices(mesh)
        verts[V] = mesh.vertices[V]
    end

    # add a vertex in each edge centroid
    # for (E, Edge) in enumerate(cells(edges))
    for (i,E) in enumerate(edges)
        verts[NV+i] = cartesian(center(chart(edges, E)))
    end

    # Build a matrix that given a coarse face gives
    # access to the coarse edges making up its boundary
    D = copy(transpose(connectivity(edges, faces, identity)))
    rows, vals = rowvals(D), nonzeros(D)

    # add four faces in each coarse face
    nf = 4NF
    fcs = zeros(celltype(mesh), nf)

    for F in 1 : numcells(faces)

        # retrieve the indices of the fine vertices in the centroids
        # of the edges in the boundary of coarse face F
        es = zeros(Int,3)
        for i in nzrange(D, F)
            E = rows[i]
            j = abs(vals[i]) # local index of edge E in face F
            @assert 1 <= j <= 3
            e = NV + E
            es[j] = e
        end

        # build the four triangles
        r1 = mesh.faces[F][1]
        r2 = mesh.faces[F][2]
        r3 = mesh.faces[F][3]

        fcs[4(F-1)+1] = index(r1,es[3],es[2])
        fcs[4(F-1)+2] = index(r2,es[1],es[3])
        fcs[4(F-1)+3] = index(r3,es[2],es[1])
        fcs[4(F-1)+4] = index(es[1],es[2],es[3])
    end

    Mesh(verts, fcs)
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
    fcs = zeros(celltype(mesh), nf)
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
