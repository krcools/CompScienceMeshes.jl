export barycentric_refinement

function barycentric_refinement{U}(mesh::Mesh{U,2})

    PointType = vertextype(mesh) #eltype(mesh.vertices)
    CellType = celltype(mesh) #eltype(mesh.faces)

    Verts = skeleton(mesh, 0)
    Edges = skeleton(mesh, 1)

    NV, NE = numcells(Verts), numcells(Edges)

    nv = NV + NE
    ne = 2NE

    verts = zeros(PointType, nv)
    for i in 1:NV
        verts[i] = mesh.vertices[i]
    end

    for i in 1:NE
        v = mesh.vertices[mesh.faces[i]]
        verts[NV+i] = (v[1] + v[2]) / 2
    end

    # build the edge to vertex lookup table
    D = connectivity(mesh, 0, Verts, Edges)
    D = transpose(D)

    edges = zeros(CellType, ne)
    rows, vals = rowvals(D), nonzeros(D)
    for E in 1:NE
        c = NV + E
        a = Edges.faces[E][1]
        b = Edges.faces[E][2]
        edges[2(E-1) + 1] = Vec(a,c)
        edges[2(E-1) + 2] = Vec(c,b)
    end

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
        v = edges.vertices[edges.faces[E]]
        verts[NV+E] = (v[1] + v[2]) / 2
    end

    # add a vertex in each face centroid
    for F in 1 : numcells(faces)
        v = faces.vertices[faces.faces[F]]
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

            fcs[6(F-1)+2(j-1)+1] = Vec(a,e,c)
            fcs[6(F-1)+2(j-1)+2] = Vec(b,c,e)
        end
    end

    Mesh(verts, fcs)
end
