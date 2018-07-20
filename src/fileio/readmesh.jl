"""
    readmesh(filename)

Reads a mesh in *in* format from `filename`. The format follows:

    1
    V C
    x1_1    x1_2    ... x1_U
    x2_1    x2_2    ... x2_U
    ...
    xV_1    xV_2    ... xV_U
    i1_1    i1_2    ... i1_D1
    i2_1    i2_2    ... i2_D1
    ...
    iC_1    iC_2    ... iC_D1

where `U` is the universedimension of the mesh, `D1` the dimension
of the mesh plus one, `V` the number of vertices, and `C` the number
of cells in the mesh.
"""
function readmesh(filename)
    open(filename) do f
        # multi-mesh files are not supported
        readline(f)

        # read the number of vertices and faces
        l = readline(f)

        sl = split(l)
        num_vertices = parse(Int, sl[1])
        num_faces    = parse(Int, sl[2])

        # Determine the universe dimension
        p = position(f)
        l = readline(f)
        udim = length(split(l))
        seek(f,p)

        T = Float64
        P = SVector{udim,T}
        vertices = zeros(P, num_vertices)
        for i = 1 : num_vertices
            l = readline(f)
            vertices[i] = P(parse.(T,split(l)))
            #vertices[i] = P(float(split(l)))
        end

        # determin the mesh dimension (plus 1)
        p = position(f)
        l = readline(f)
        dim1 = length(split(l))
        seek(f,p)

        # TODO: remove explicit reference to dimension
        #I = SVector{3,Int}
        C = SVector{dim1,Int}
        faces = zeros(C, num_faces)
        for i in 1 : num_faces
            l = readline(f)
            faces[i] = C([parse(Int,s) for s in split(l)])
        end

        Mesh(vertices, faces)
    end
end



"""
    writemesh(mesh, filename)

Write `mesh` to `filename` in the *in* format (see `readmesh`).
"""
function writemesh(mesh, filename)
    dl = '\t'
    open(filename, "w") do f
        println(f, 1)
        println(f, numvertices(mesh), dl, numcells(mesh))
        for v in mesh.vertices
            println(f, v[1], dl, v[2], dl, v[3])
        end
        for c in mesh.faces
            println(f, c[1], dl, c[2], dl, c[3])
        end
    end
end
