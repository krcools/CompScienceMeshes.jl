"""
    read_gmsh_mesh(filename) -> mesh::Mesh
"""
function read_gmsh_mesh(fn::AbstractString; kwargs...)
    open(fn, "r") do io
        read_gmsh_mesh(io; kwargs...)
    end
end

"""
    read_gmsh_mesh(iostream) -> mesh::Mesh

Reads the mesh nodes and elements stored in the input .msh file (`io`, output by GMSH)
    into arrays of node vectors and vertex vectors respectively.

Returns an object `mesh` of type `Mesh`, comprising both vector arrays.
"""
function read_gmsh_mesh(io; physical=nothing, dimension=2, sort=true, T=Float64)

    entity_tag = 0
    if physical != nothing
        @assert physical isa String
        while true
            thisLine = io |> readline |> strip
            thisLine == "\$PhysicalNames" && break
            thisLine == "\$Nodes" && error("Mesh file does not define physical entities")
        end
        thisLine = io |> readline |> strip
        s = split(thisLine)
        NP = parse(Int,s[1])
        print("Mesh file defined $(NP) physical entities: ")
        for i in 1:NP
            thisLine = io |> readline |> strip
            s = split(thisLine)
            entity_name = s[3]
            # @show entity_name
            if entity_name == "\"$(physical)\""
                entity_tag = parse(Int, s[2])
                println("Target entity \"$(physical)\" is $i out of $(NP).")
                break
            end
        end
        entity_tag == 0 && error("Specified physical entitiy not found")
    end

    thisLine = io |> readline |> strip
    while thisLine != "\$Nodes"
        thisLine = io |> readline |> strip
    end

    thisLine = io |> readline |> strip
    s = split(thisLine)
    NV = parse(Int, s[1])

    P = SVector{3,T}
    v = Vector{P}(undef,NV)
    for i in 1:NV
        thisLine = io |> readline |>  strip
        d = readdlm(IOBuffer(thisLine), T)
        v[i] = P(d[2], d[3], d[4])
    end

    while thisLine != "\$Elements"
        thisLine = io |> readline |> strip
    end

    thisLine = io |> readline |> strip
    s = split(thisLine)
    NF = parse(Int, s[1])

    dim2type = [1,2,4]
    type = dim2type[dimension]

    thisLine = io |> readline |> strip
    I = SVector{dimension+1,Int}
    f = Vector{I}(undef,NF)
    i = 0
    while thisLine != "\$EndElements"
        d = readdlm(IOBuffer(thisLine), Int)
        thisLine = io |> readline |> strip
        d[2] == type || continue
        entity_tag == d[4] || entity_tag == 0 || continue
        # f[i +=1] = SVector(d[end-2], d[end-1], d[end-0])
        f[i+=1] = d[end-dimension:end]
    end
    resize!(f,i)

    # @show length(v)
    # @show length(f)
    
    if sort
        # @info "sorting..."
        f = sort_sfc(v,f)
    end
    # Q = Pt{3,Float64}
    # ctrs = Q[]
    # for tr in f
    #         ctr = sum(v[tr])/3
    #         push!(ctrs,ctr)
    # end
    # sorted = sort_sfc(ctrs)
    # f = f[sorted]



    return Mesh(v, f)

end


function sort_sfc(v, f)
    T = eltype(eltype(v))
    Q = Pt{3,T}
    ctrs = Q[]
    for tr in f
            ctr = sum(v[tr])/3
            push!(ctrs,ctr)
    end
    sorted = sort_sfc(ctrs)
    return f[sorted]
end