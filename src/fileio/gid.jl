"""
    load_gid_mesh(filename) -> mesh
"""
function load_gid_mesh(fn::AbstractString,T=Float64)
    open(fn, "r") do io
        load_gid_mesh(io,T)
    end
end

"""
    load_gid_mesh(stream) ->mesh
"""
function load_gid_mesh(file,T)

    lineCount = countlines(file);
    seekstart(file);
    readline(file); # MESH dimension 3 ElemType Triangle Nnode 3
    readline(file); # Coordinates

    # Read the vertices
    vertices = Vector{Pt{3,T}}(undef,lineCount)
    vertexCount = 0;
    while !eof(file)
            line = split(readline(file));
            size(line, 1) == 4 || break; # End Coordinates
            vertexCount += 1;
            #vertices[vertexCount] = Pt{3,Float64}(map(Meta.parse, line[2:4]));
            vertices[vertexCount] = Pt{3,T}(parse.(T,line[2:4]))
    end
    vertices= vertices[1:vertexCount]

    readline(file); # "\n"
    readline(file); # Elements

    # Read the triangles
    triangles = Vector{Pt{3,Int}}(undef,lineCount-vertexCount);
    triangleCount = 0;
    while !eof(file)
            line = split(readline(file));
            size(line, 1) == 4 || break; # End Elements
            triangleCount += 1;
            triangles[triangleCount] = Pt{3,Int}(map(Meta.parse, line[2:4]));
    end
    triangles = triangles[1:triangleCount]
    # Mesh(vertices[1:vertexCount], triangles[1:triangleCount])

    # Place them on a space filling curve
    Q = Pt{3,T}
    ctrs = Q[]
    for tr in triangles
            ctr = sum(vertices[tr])/3
            push!(ctrs,ctr)
    end
    sorted = sort_sfc(ctrs)
    triangles = triangles[sorted]

    Mesh(vertices, triangles)
end
