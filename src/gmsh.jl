"""
    read_gmsh_mesh(filename) -> mesh::Mesh
"""
function read_gmsh_mesh(fn::AbstractString)
    open(fn, "r") do io
        read_gmsh_mesh(io)
    end
end

"""
    read_gmsh_mesh(iostream) -> mesh::Mesh
"""
function read_gmsh_mesh(io)

    thisLine = io |> readline |> strip
    while thisLine != "\$Nodes"
        thisLine = io |> readline |> strip
    end

    thisLine = io |> readline |> strip
    s = split(thisLine)
    NV = parse(Int, s[1])

    P = Vec{3,Float64}
    v = Vector{P}(NV)
    for i in 1:NV
        thisLine = io |> readline |>  strip
        d = readdlm(IOBuffer(thisLine), Float64)
        v[i] = P(d[2], d[3], d[4])
    end

    while thisLine != "\$Elements"
        thisLine = io |> readline |> strip
    end

    thisLine = io |> readline |> strip
    s = split(thisLine)
    NF = parse(Int, s[1])

    thisLine = io |> readline |> strip
    f = Vector{Vec{3,Int}}(NF); i = 0
    while thisLine != "\$EndElements"
        d = readdlm(IOBuffer(thisLine), Int)
        thisLine = io |> readline |> strip
        d[2] == 2 || continue
        f[i +=1] = Vec(d[end-2], d[end-1], d[end-0])
    end
    resize!(f,i)

    return Mesh(v, f)

end
