function load_gid_mesh(fn::AbstractString)
    open(fn, "r") do io
        load_gid_mesh(io)
    end
end

function load_gid_mesh(file)
    #file = open(fileName, "r");

    lineCount = countlines(file);
    seekstart(file);
    readline(file); # MESH dimension 3 ElemType Triangle Nnode 3
    readline(file); # Coordinates

    # Read the vertices
    vertices = Vector{FixedSizeArrays.Vec{3,Float64}}(lineCount);
    vertexCount = 0;
    while !eof(file)
            line = split(readline(file));
            size(line, 1) == 4 || break; # End Coordinates
            vertexCount += 1;
            vertices[vertexCount] = FixedSizeArrays.Vec{3,Float64}(map(parse, line[2:4]));
    end

    readline(file); # "\n"
    readline(file); # Elements

    # Read the triangles
    triangles = Vector{FixedSizeArrays.Vec{3,Int}}(lineCount-vertexCount);
    triangleCount = 0;
    while !eof(file)
            line = split(readline(file));
            size(line, 1) == 4 || break; # End Elements
            triangleCount += 1;
            triangles[triangleCount] = FixedSizeArrays.Vec{3,Int}(map(parse, line[2:4]));
    end

    #close(file);

    return Mesh{3,3, Float64}(vertices[1:vertexCount], triangles[1:triangleCount])
end
