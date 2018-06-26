"""
    Loop_subdivsion(mesh) -> refinement

Construct a refinement of `mesh` by Loop subdivision scheme.
Every face is subdived into four small faces and use weights to smooth the surface.

Only defined for 2D meshes.
"""

function Loop_subdivision(mesh::Mesh{U,3}) where U

    D1 = 3

    edges = skeleton(mesh,1)
    faces = skeleton(mesh,2)

    connectMat = connectivity(edges,faces)

    NV, NE, NF = numvertices(mesh), numcells(edges), numcells(faces)

    nv = NV + NE

    verts = Array{vertextype(mesh)}(nv)

    # copy over the existing vertices
    for V in 1 : numvertices(mesh)
        # create a mask for new vertex point
        RingverticesIndices = []
        # loop over all the faces
        for (F, Face) in enumerate(cells(faces))
            for iv = 1 : 3
                # if V is a vertex of the face add all the vertices to the mask
                if Face[iv] == V
                    append!(RingverticesIndices, Face)
                end
            end
        end
       # delete duplicate vertices
        RingverticesIndices = unique(RingverticesIndices)
       # calculate valence for vertex V
        valence = length(RingverticesIndices) - 1
        # calculate weights
        if valence != 3
            weight = 3.0 / (8.0 * valence)
        else
            weight = 3.0 / 16.0
        end
        newvert = [0.0, 0.0, 0.0]
        for iv = 1 : length(RingverticesIndices)
             if RingverticesIndices[iv] == V
                     newvert += mesh.vertices[RingverticesIndices[iv]] * (1.0 - valence * weight)
             else
                     newvert += mesh.vertices[RingverticesIndices[iv]] * weight
             end
        end
        verts[V] = newvert
    end

    # add a vertex for each edge
    for (E, Edge) in enumerate(cells(edges))
        # create a mask for new edge point
        EndVerticesIndices = Edge
        FacesIndex = []
        RelatedVerticesIndices = []
        # Loop over all faces
        for (F, Face) in enumerate(cells(faces))
            if connectMat[F,E] != 0
                # Find out two faces in the mask
                append!(FacesIndex, F)
                ## Loop over vertices in the faces find out related vertices
                for iv = 1 : length(Face)
                    if Face[iv] != Edge[1] && Face[iv] != Edge[2]
                            append!(RelatedVerticesIndices, Face[iv])
                    end
                end
            end
        end
        # verts[NV+E] = cartesian(center(chart(edges, Edge)))
        # compute new edge vertices
        if length(RelatedVerticesIndices) != 2 || length(EndVerticesIndices) != 2
            print("The number of vertices in the mask for a edge point is not right!")
        end
        newvert = [0.0, 0.0, 0.0]
        for iv = 1 : length(EndVerticesIndices)
            newvert = newvert + mesh.vertices[EndVerticesIndices[iv]] * 0.375
        end
        for iv = 1: length(RelatedVerticesIndices)
            newvert = newvert + mesh.vertices[RelatedVerticesIndices[iv]] * 0.125
        end
        verts[NV+E] = newvert
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
