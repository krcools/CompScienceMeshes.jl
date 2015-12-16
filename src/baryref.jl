export barycentric_refinement

function barycentric_refinement{T}(mesh::Mesh{3,3,T})

    U, D1 = 3, 3

    edges = cells(mesh,1)
    faces = cells(mesh,2)

    NV, NE, NF = numvertices(mesh), numcells(edges), numcells(faces)

    nv = NV + NE + NF
    verts = Array(Point{U,T}, nv)
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

    D = transpose(connectivity(mesh, 1, edges, faces, op=identity))
    rows, vals = rowvals(D), nonzeros(D)

    # add six faces in each coarse face
    nf = 6NF
    fcs = zeros(Vec{D1,Int}, nf)
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
