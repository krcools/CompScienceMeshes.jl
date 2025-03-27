using StaticArrays

"""
    meshicosphere(;nu::Int = 1, r::F = 1.0, nr_verts) where F

returns Mesh(vertices, faces)

Function gives a geodesic icosahedron with subdivision frequency nu. Frequency
nu = 1 returns regular unit icosahedron, and nu>1 preformes subdivision. 
The result is a structured sphere.

Using subdivision frequency, possible number of vertices grows with nu as
    [12+10*(nu+1)*(nu-1) for nu in range(1,33)]
    which gives 
    [12, 42, 92, 162, 252, 362, 492, 642, 812, 1002, 1212, 1442, 1692, 1962, 
     2252, 2562, 2892, 3242, 3612, 4002, 4412, 4842, 5292, 5762, 6252, 6762, 
     7292, 7842, 8412, 9002, 9612, 10242]

For example, nu = 3, divides each triangle of the polygon of icosahedron to

                    .
                   / \\
                  .---.
                 /\\ / \\
                .---.---.
               /\\ /\\ / \\
              .---.---.---.

kwargs:
Radius - is scaled with respect to the unit icosphere from given nu. 
nr_verts -
If nr_verts is given, nu will be adjusted such that icosphere contains
at least nr_verts vertices.  

Original author and python code by vand@dtu.dk, 2014, 2017, 2021. Developed 
in connection with https://ieeexplore.ieee.org/document/7182720

"""
function meshicosphere(nu::Int = 1, r::F = 1.0, nr_verts = nothing) where F
    #Unit icosahedron for radius r = 1 (icosahedron with no subdivisions)
    vertices, faces = icosahedron(r)
    if !isnothing(nr_verts)
        nu_min = Int(ceil(sqrt(max(1 + (nr_verts - 12)/10, 1))))
        nu = max(nu, nu_min)
    end

    if nu > 1
        vertices, faces = subdivide_mesh(vertices, faces, nu)
        for i in 1 : length(vertices)
            vertices[i] = r.*vertices[i]./sqrt.(sum(x->x.^2, vertices[i]))
        end
    end 
    return Mesh(vertices, faces)
end

"""
    icosahedron(r::F) where F

Returns a unit centred icosahedron
"""
function icosahedron(r::F) where F

    ϕ = (1 + sqrt(5))/2
    vertices = Vector{SVector{3, F}}([
        [0, 1, ϕ], [0, -1, ϕ], [1, ϕ, 0], [-1, ϕ, 0], 
        [ϕ, 0, 1], [-ϕ, 0, 1]
        ]*(r/sqrt(1 + ϕ^2))
        )
    vertices = vcat(vertices, -vertices)
    faces = Vector{SVector{3, Int}}([
        [1, 6, 2], [1, 4, 6], [1, 3, 4], [1, 5, 3], [1, 2, 5],
        [2, 6, 9], [6, 4, 11], [4, 3, 8], [3, 5, 12], [5, 2, 10],
        [8, 12, 7], [12, 10, 7], [10, 9, 7], [9, 11, 7], [11, 8, 7],
        [3, 12, 8], [5, 10, 12], [2, 9, 10], [6, 11, 9], [4, 8, 11]
        ])
    return vertices, faces
end

"""
    subdivide_mesh(vertices, faces, nu)

Function divides the unit isocahedron along each triangle
nu denotes the number of elements each side of triangle is divided
into, in a icosahedron.
"""
function subdivide_mesh(vertices, faces, nu)
    edges = []
    #edges collects edges AB, AC, BC for all faces
    for i in 1 : length(faces)
       append!(edges, [
            sort([faces[i][1], faces[i][2]]), 
            sort([faces[i][2], faces[i][3]]), 
            sort([faces[i][1], faces[i][3]])
            ])
    end
    #no common edge is repeated
    unique!(sort!(edges))

    F = length(faces)
    V = length(vertices)
    E = length(edges)

    subfaces = Vector{SVector{3, Int}}(undef, F*nu^2)
    #= total number of (all) vertices === subvertices would be 
    the number of counted vertices, plus
    the number of (sub-)vertices on each edge, plus the 
    number of vertices on the inner faces =#
    T = typeof(vertices[1])
    subvertices = Vector{T}(
        undef, 
        Int(V + E*(nu - 1) + F*floor((nu - 1)*(nu - 2)/2))
        )
    subvertices[begin:V] = vertices #vertices of the unit icosahedron

    edge_indices = Dict()
    #assigning keys to the edges
    for i in 1 : V
        edge_indices[i] = Dict()
    end
    for i in 1 : E
        edge_indices[edges[i][1]][edges[i][2]] = i - 1
        edge_indices[edges[i][2]][edges[i][1]] = -i + 1
    end

    template = faces_template(nu)
    ordering = vertex_ordering(nu)

    reordered_template = Vector{Int}[]
    for i in eachindex(template)
        append!(reordered_template, [ordering[template[i]]])
    end
    #w gives the weights of the weights of the vertices along the edge
    w = collect((1 : nu - 1)./nu)
    for e in 1 : E #e is the number of edge, 
        #(e - 1) would count the number of counted edges
        edge = edges[e]
        for k in 1 : nu - 1  #k gives the number of subvertices on the edge
            subvertices[V + (e - 1)*(nu - 1) + k] = (
                w[end - k + 1]*vertices[edge[1]] 
                + w[k]*vertices[edge[2]]
                )
        end
    end

    r = collect(1 : nu - 1)
    for f in 1 : F
        T = collect(range(Int(floor(
            (f - 1)*(nu - 1)*(nu - 2)/2)) + E*(nu - 1) + V + 1, 
            Int(floor((f)*(nu - 1)*(nu - 2)/2)) + E*(nu - 1) + V, step=1)
            )
        eAB = edge_indices[faces[f][1]][faces[f][2]]
        eAC = edge_indices[faces[f][1]][faces[f][3]]
        eBC = edge_indices[faces[f][2]][faces[f][3]]
        AB = rreverse((abs(eAB))*(nu - 1) + (V) .+ r, eAB < 0)
        AC = rreverse((abs(eAC))*(nu - 1) + (V) .+ r, eAC < 0)
        BC = rreverse((abs(eBC))*(nu - 1) + (V) .+ r, eBC < 0)
        VEF = vcat(faces[f], AB, AC, BC, T)
        for i in 1:nu^2
            subfaces[(f - 1)*nu^2 + i] = [
                VEF[reordered_template[i][1]], 
                VEF[reordered_template[i][2]], 
                VEF[reordered_template[i][3]]
                ]
        end
        subvertices[T] = inside_points(subvertices[AB], subvertices[AC])
    end
    return (subvertices, subfaces)
end

"""
    rreverse(vector, flag)

Function reverses the order of the vector
A flag is raised when the order of the vertices along an edge
is not anticlockwise or the edge direction is negative.
"""
function rreverse(vector, flag)
    if length(vector) == 1
        skip
    elseif flag == true
        vector = Base.reverse(vector)
    end
    return vector
end

function faces_template(nu)
    faces = Vector{Int}[]
    for i in 1 : nu
        vertex = Int(floor(i*(i - 1)/2))
        skip = i
        for j in 1 : i - 1
            append!(faces, [
                [j + vertex, j + vertex + skip, j + vertex + skip + 1]
                ])
            append!(faces, [
                [j + vertex, j + vertex + skip + 1, j + vertex + 1]
                ])
        end
        #the last face that is unpaired
        append!(faces, [
            [i + vertex, i + vertex + skip, i + vertex + skip + 1]
            ])
    end
    return faces
end

function vertex_ordering(nu)
    left = [j for j in 4 : nu + 2]
    right = [j for j in nu + 3 : 2*nu + 1]
    bottom = [j for j in 2*nu + 2 : 3*nu]
    inside = [j for j in 3*nu + 1 : Int(floor((nu + 1)*(nu + 2)/2))]
    o = [1] #topmost corner
    for i in 1 : nu - 1
        append!(o, left[i])
        append!(o, inside[
            Int(floor((i - 1)*(i - 2)/2)) + 1 : Int(floor((i - 1)*i/2))
            ])
        append!(o, right[i])
    end
    append!(o, 2)
    append!(o, bottom)
    append!(o, 3)
    return o
end

function inside_points(vAB, vAC)
    l = length(vAB)
    v = []
    if l == 1
        skip
    else
        for i in 1 : length(vAB) - 1
            w = collect((1 : i)./(i + 1))     
            for k in 1 : i
                append!(v, [
                    SVector(w[end - k + 1].*vAB[i + 1][:] + w[k].*vAC[i + 1][:])
                    ])
            end
        end
    end
    return v
end