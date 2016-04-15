function read_gmsh_mesh(io)

    thisLine = io |> readline |> strip
    while thisLine != "\$Nodes"
        thisLine = io |> readline |> strip
    end

    thisLine = io |> readline |> strip
    s = split(thisLine)
    NV = parse(Int, s[1])

    #v = readdlm(io, Float64; dims=(NV,4))
    v = Vector{Point{3,Float64}}(NV)
    for i in 1:NV
        thisLine = io |> readline |>  strip
        d = readdlm(IOBuffer(thisLine), Float64)
        v[i] = Point(d[2], d[3], d[4])
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
