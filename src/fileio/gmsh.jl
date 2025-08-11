"""
    load_gmsh_mesh(meshfile; element=:triangle, vertextype=Float64, order=1, sort=true)

Reads the mesh nodes and elements stored in the input .msh file (`io`, output by Gmsh) to generated
a mesh.
If `udim=2`, then the z-component is removed resulting a mesh for 2D BEM simulations
(we cannot deduce this from element=:line, because these could be used for thin wire simulations)
Only elements of the specified type `element` and `order` are loaded.
When `sort==true` (default), the vertices are sorted such that the lie on a Hilbert space filling
curve.

TODO: Curvilinear meshes are stored, for the time being, in the same mesh type.
"""
function load_gmsh_mesh(meshfile;
    element=:triangle,
    udim=3,
    vertextype=Float64,
    order=1,
    physical=nothing,
    sort=true
)

    gmsh.initialize()

    try
        gmsh.open(meshfile)
    catch e
        gmsh.finalize()
        error("Mesh file could not be opened!")
    end
    # Get all nodes
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes(
        -1, # dim<0: Get nodes of entities of any dimension (all node)
        -1, # tag<0: Get the nods of all entities of dimension dim
        true, # includeBoundary (default false)
        true #  returnParametricCoord (default true)
    )

    if physical===nothing
        element_types, element_tags, element_node_tags = gmsh.model.mesh.getElements()
    else
        # Get tuples of the form (dim, tag)
        # dim: topological dimension (0: point, 1: curve, 2: surfaces, 3: volumes)
        # tag: the unique tag number
        # For example: 1 2 "Coast" --> dim = 1, tag = 2
        phys_groups = gmsh.model.getPhysicalGroups()

        # Each item is (dim, tag)
        phys_name_to_tag = Dict(gmsh.model.getPhysicalName(dim, tag) => (dim, tag)
            for (dim, tag) in phys_groups)

        @show (dim, phys_tag) = phys_name_to_tag[physical]

        # Get the geometrical entities (e.g., surfaces or volumes) belonging to that physical group
        @show entities = gmsh.model.getEntitiesForPhysicalGroup(dim, phys_tag)

        # Get elements for these entities
        element_types = Int[]
        element_tags = Vector{Int}[]
        element_node_tags = Vector{Int}[]

        for ent in entities
            elemtypes, elemtags, nodetags = gmsh.model.mesh.getElements(dim, ent)
            append!(element_types, elemtypes)
            append!(element_tags, vcat(elemtags))
            append!(element_node_tags, vcat(nodetags))
        end

    end

    gmsh.finalize()

    # While I added a feature that this routine trashes automatically the z-component,
    # I am not using it here in order to be able to use the sort function
    vertices = gmsh_to_mesh_nodes(node_tags, node_coords; udim=3, vertextype=vertextype)
    elements = gmsh_get_elements(element_types, element_tags, element_node_tags; element=element, order=order)

    if sort
        # @info "Sorting vertices to lie on a Hilbert space filling curve"
        elements = sort_sfc(vertices, elements)
    end

    if udim == 2
        vertices2d = [SVector(v[1], v[2]) for v in vertices]
        vertices = vertices2d
    end

    # TODO: once we support hexahedrons, we need to catch that one too
    if element == :quadrangle
        return QuadMesh(vertices, elements)
    else
        return Mesh(vertices, elements)
    end
end

"""
    gmsh_to_mesh_nodes(node_tags, node_coords; udim=3, vertextype=Float64)

Given the vectors node_tags and node_coords provided by Gmsh's getNodes() function,
return a vector containing with the nodes' coordinates.
Use the keyword argument `vertextype=Float32` to reduce the precision and use `udim=2`
to ensure that the mesh is loaded a 2D-mesh by trashing the z-coordinate.
"""
function gmsh_to_mesh_nodes(node_tags, node_coords; udim=3, vertextype=Float64)
    maxnodeindex = maximum(node_tags)

    SV = SVector{udim,vertextype}

    if udim==2 || udim==3
        NonInitVertex = SVector{udim, vertextype}(fill(NaN, udim)...)
    else
        error("Unsupported universe dimension")
    end

    # It is unclear to me if there is a chance
    # the there are node numbers appearing larger
    # then the number of nodes because Gmsh skips some
    # indices
    # In this case, the skipped elements are set to NaN-vectors
    # ALTERNATIVELY: We could not add those nodes at all
    # but then we need to correct all node indices
    vertices = fill(NonInitVertex, maxnodeindex) #SVector#Vector{SV}(SVe, maxnodeindex)

    if length(node_coords) == 3*length(node_tags)
        dimfactor = 3
    elseif length(node_coords) == 2*length(node_tags)
        dimfactor = 2
    else
        error("The coordinate vector does not match the number of vertices")
    end

    for (i, node) in enumerate(node_tags)
        x = vertextype(node_coords[dimfactor*(i-1)+1])
        y = vertextype(node_coords[dimfactor*(i-1)+2])

        if udim==2
            coord = SV(x, y)
        elseif udim==3
            z = vertextype(node_coords[3(i-1)+3])
            coord = SV(x, y, z)
        else
            error("Unsupported universe dimension")
        end

        vertices[node] = coord
    end

    return vertices
end

"""
    gmsh_get_elements(element_types, element_tags, element_node_tags; element=:triangle, order=1)

Given the vectors element_types, element_tags, and element_node_tags
provided by Gmsh's getElements() function, return a vector containing the elements
of type `element` and `order`.
An element (line, triangle, hexahedron etc.) is described by an SVector containing
the node indices of the interpolation nodes. We use the same ordering as provided
by Gmsh.
"""
function gmsh_get_elements(element_types, element_tags, element_node_tags; element=:triangle, order=1)

    # The Gmsh Element Type ID (second argument of the tuple) is from
    # https://gitlab.onelab.info/gmsh/gmsh/blob/master/src/common/GmshDefines.h
    # Only some of the element IDs are explained in the manual
    # Format (Number of interpolation points, Gmsh ID)
    line = SVector(
        (2, 1),     # Order 1: 2-node line
        (3, 8),     # Order 2: 3-node line
        (4, 26),    # Order 3: 4-node line
        (5, 27),    # Order 4: 5-node line
        (6, 28),    # Order 5: 6-node line
        (7, 62),    # Order 6: 7-node line
        (8, 63),    # Order 7: 8-node line
        (9, 64),    # Order 8: 9-node line
        (10, 65),   # Order 9: 10-node line
        (11, 66)    # Order 10: 11-node line
    )

    # Full order
    triangle = SVector(
        (3, 2),     # Order 1: 3-node triangle
        (6, 9),     # Order 2: 6-node triangle
        #(9, 20),
        (10, 21),   # Order 3: 10-node triangle
        #(12, 22),
        (15, 23),   # Order 4: 15-node triangle
        (21, 25),   # Order 5: 21-node triangle
        (28, 42),   # Order 6: 28-node triangle
        (36, 43),   # Order 7: 36-node triangle
        (45, 44),   # Order 8: 45-node triangle
        (55, 45),   # Order 9: 55-node triangle
        (66, 46)    # Order 10: 66-node triangle
        #(18, 52),
        #(24, 54),
        #(27, 55),
        #(30, 56)
    )

    # Full order (no serendipity elements)
    quadrangle = SVector(
        (4, 3),     # Order 1: 4-node quadrangle
        (9, 10),    # Order 2: 9-node quadrangle
        (16, 36),   # Order 3: 16-node quadrangle
        (25, 37),   # Order 4: 25-node quadrangle
        (36, 38),   # Order 5: 36-node quadrangle
        (49, 47),   # Order 6: 49-node quadrangle
        (64, 48),   # Order 7: 64-node quadrangle
        (81, 49),   # Order 8: 81-node quadrangle
        (100, 50),  # Order 9: 100-node quadrangle
        (121, 51)   # Order 10: 121-node quadrangle
    )

    # Full order
    tetrahedron = SVector(
        (4, 4),       # Order 1: 4-node tetrahedron
        (10, 11),     # Order 2: 10-node tetrahedron
        (20, 29),     # Order 3: 20-node tetrahedron
        (35, 30),     # Order 4: 35-node tetrahedron
        (56, 31),     # Order 5: 56-node tetrahedron
        (84, 71),     # Order 6: 84-node tetrahedron
        (120, 72),    # Order 7: 120-node tetrahedron
        (165, 73),    # Order 8: 165-node tetrahedron
        (220, 74),    # Order 9: 220-node tetrahedron
        (286, 75)     # Order 10: 286-node tetrahedron
    )

    # Full order
    hexahedron = SVector(
        (8, 5),        # Order 1: 8-node hexahedron
        (27, 12),      # Order 2: 27-node hexahedron
        (64, 28),      # Order 3: 64-node hexahedron
        (125, 30),     # Order 4: 125-node hexahedron
        (216, 32),     # Order 5: 216-node hexahedron
        (343, 34),     # Order 6: 343-node hexahedron
        (512, 75),     # Order 7: 512-node hexahedron
        (729, 76),     # Order 8: 729-node hexahedron
        (1000, 77),    # Order 9: 1000-node hexahedron
        (1331, 78)     # Order 10: 1331-node hexahedron
    )

    if element == :line
        gmsh_element_id = line[order][2]
        numnodes = line[order][1]
    elseif element == :triangle
        gmsh_element_id = triangle[order][2]
        numnodes = triangle[order][1]
    elseif element == :quadrangle
        gmsh_element_id = quadrangle[order][2]
        numnodes = quadrangle[order][1]
    elseif element == :tetrahedron
        gmsh_element_id = tetrahedron[order][2]
        numnodes = tetrahedron[order][1]
    elseif element == :hexahedron
        gmsh_element_id = hexahedron[order][2]
        numnodes = hexahedron[order][1]
    else
        error("$element is not supported")
    end

    position_in_gmsh_arrays = -1
    for (i, type) in enumerate(element_types)
        if type == gmsh_element_id
            position_in_gmsh_arrays = i
            break
        end
    end

    SV = SVector{numnodes, Int64}
    NonInitElement = SV(fill(0, numnodes)...) # if its norm is zero, we know it is a none-entry

    if position_in_gmsh_arrays == -1
        # If == -1, then we know that this element type / order combination does not exist
        return SV[] # Return empty array
    end

    # TODO: Currently, we are not using the element tags
    # This might be needed if not only type of element is loaded
    # but, for example, tetrahedrons and triangles.
    elementtags = element_tags[position_in_gmsh_arrays]
    #maxelementindex = maximum(elementtags)

    element_nodes = element_node_tags[position_in_gmsh_arrays]
    elements = fill(NonInitElement, length(elementtags))

    for i in eachindex(elementtags)
        elements[i] = SV(element_nodes[numnodes*(i-1)+1:numnodes*i]...)
    end

    return elements
end


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