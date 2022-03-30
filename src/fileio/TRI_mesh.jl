"""
    read_TRI_mesh(filename) -> mesh::Mesh
"""
function read_TRI_mesh(fn::AbstractString)
    open(fn, "r") do io
        read_TRI_mesh(io)
    end
end

"""
    read_TRI_mesh(mesh_filename) -> mesh::Mesh

Imports a surface mesh (stored in an ASCII file named `mesh_filename`) into a
    `Mesh` object (i.e. node list and element list).

NOTE: The contents of `mesh_filename` must include the file extension, and the
    file must be stored in the current directory.

"""
function read_TRI_mesh(io,T=Type::Float64)

    thisLine = io |> readline |> strip
    while thisLine != "POINTS"
        thisLine = io |> readline |> strip
    end

    thisLine = io |> readline |> strip
    s = split(thisLine)
    NV = parse(Int, s[1])

    P = SVector{3,T}
    v = Vector{P}(NV)
    for i in 1:NV
        thisLine = io |> readline |>  strip
        d = readdlm(IOBuffer(thisLine), T)
        v[i] = P(d[2], d[3], d[4])
    end

    while thisLine != "TRIANGLES"
        thisLine = io |> readline |> strip
    end

    thisLine = io |> readline |> strip
    s = split(thisLine)
    NF = parse(Int, s[1])

    thisLine = io |> readline |> strip
    P = SVector{3,Int64}
    f = Vector{P}(NF); i = 0
    while thisLine != "END_TRIANGLES"
        d = readdlm(IOBuffer(thisLine), Int64)
        thisLine = io |> readline |>  strip
        f[i += 1] = P(d[2], d[3], d[4])
    end

    return Mesh(v, f)

end
