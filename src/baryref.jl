

"""
    barycentric refinement(mesh) -> refined_mesh

Create the mesh obtained by inserting an extra vertex in the barycenters of all cells and
recusively creating fine cells by connecting the barycenter of a k-cell to the already
constructed refined (k-1)-cells on its boundary.
"""
function barycentric_refinement(mesh::Mesh{U,2}) where U

    # Get the points and the faces (segments) we will use in the refinement process
    Edges = skeleton(mesh, 1)

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
    for (i, Face) in enumerate(cells(mesh))
        verts[NV+i] = cartesian(center(chart(mesh, Face)))
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

function barycentric_refinement(mesh::Mesh{U,3}) where U

    D1 = 3

    edges = skeleton(mesh,1)
    faces = skeleton(mesh,2)

    NV, NE, NF = numvertices(mesh), numcells(edges), numcells(faces)

    nv = NV + NE + NF
    verts = Array{vertextype(mesh)}(undef,nv)
    for V in 1 : numvertices(mesh)
        verts[V] = mesh.vertices[V]
    end

    # add a vertex in each edge centroid
    for (E,Edge) in enumerate(cells(edges))
        verts[NV+E] = cartesian(center(chart(edges, Edge)))
    end

    # add a vertex in each face centroid
    for (F,Face) in enumerate(cells(faces))
        verts[NV+NE+F] = cartesian(center(chart(faces, Face)))
    end

    D = copy(transpose(connectivity(edges, faces, identity)))
    rows, vals = rowvals(D), nonzeros(D)

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

    Mesh(verts, fcs)
end


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
    for (E, Edge) in enumerate(cells(edges))
        verts[NV+E] = cartesian(center(chart(edges, Edge)))
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
