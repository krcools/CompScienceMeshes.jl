export barycentric_refinement

function barycentric_refinement{U}(mesh::Mesh{U,2})
    # We will use the following types to create the longer array of vertices and faces after the refienment
    PointType = vertextype(mesh) #eltype(mesh.vertices)
    CellType = celltype(mesh) #eltype(mesh.faces)
    # Get the points and the faces (segments) we will use in the refinement process
    Verts = skeleton(mesh, 0)  # we didn't use this though
    Edges = skeleton(mesh, 1)
    # Note thier number and use it to find the new number of points and segments(faces) after refienmnts
    NV, NE = numcells(Verts), numcells(Edges)

    nv = NV + NE
    ne = 2NE
    # define the new vertices array that will hold the old points(coarse) and the new points(finer mesh) of vertices
    verts = zeros(PointType, nv)
    # Now we assign the coarse mesh vetrices (old points) in the upper half excatly as the original coarse mesh
    for i in 1:NV
        verts[i] = mesh.vertices[i]
    end
    # for the second half we store the refienment points and we calculate them using the old points
    # We insert a point between every pair of points in the original coarse mesh by (point1+ point2)\2
    for i in 1:NE
        #v = mesh.vertices[mesh.faces[i]]
        v = cellvertices(mesh,i)
        verts[NV+i] = (v[1] + v[2]) / 2
    end
    # Now we create the faces using the new points as well
    #example
      # oringinal points_index (1,2,3,4) , original faces ((1,2), (2,3), (3,4),(4,5))
      # new points_index(1,2,3,4,5,6,7,8), new_faces((1,5),(5,2),(2,6),(6,3),(3,7),(7,4),(4,8),(8,1))
      #                  a,b      For       E=1 ->   (a,c),(c,b) and so on
    edges = zeros(CellType, ne)

    for E in 1:NE
        c = NV + E
        a = Edges.faces[E][1]
        b = Edges.faces[E][2]
        edges[2(E-1) + 1] = Vec(a,c)
        edges[2(E-1) + 2] = Vec(c,b)
    end
    # return the new mesh after refienments
    return Mesh(verts, edges)
end

function barycentric_refinement{U}(mesh::Mesh{U,3})

    D1 = 3

    edges = skeleton(mesh,1)
    faces = skeleton(mesh,2)

    NV, NE, NF = numvertices(mesh), numcells(edges), numcells(faces)

    nv = NV + NE + NF
    verts = Array{vertextype(mesh)}(nv)
    for V in 1 : numvertices(mesh)
        verts[V] = mesh.vertices[V]
    end

    # add a vertex in each edge centroid
    for E in 1 : numcells(edges)
        #v = edges.vertices[edges.faces[E]]
        v = cellvertices(edges, E)
        verts[NV+E] = (v[1] + v[2]) / 2
    end

    # add a vertex in each face centroid
    for F in 1 : numcells(faces)
        #v = faces.vertices[faces.faces[F]]
        v = cellvertices(faces, F)
        verts[NV+NE+F] = (v[1] + v[2] + v[3]) / 3
    end

    #D = transpose(connectivity(mesh, 1, edges, faces, op=identity))
    D = transpose(connectivity(edges, faces, identity))
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

            # fcs[6(F-1)+2(j-1)+1] = Vec(a,e,c)
            # fcs[6(F-1)+2(j-1)+2] = Vec(b,c,e)
            fcs[6(F-1)+2(j-1)+1] = index(a,e,c)
            fcs[6(F-1)+2(j-1)+2] = index(b,c,e)
        end
    end

    Mesh(verts, fcs)
end
