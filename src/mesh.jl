export numcells
export numvertices, dimension, vertices
export meshsegment, meshrectangle, meshcircle
export readmesh, boundary, vertextocellmap, cells, connectivity, buildvolmpairs
export relorientation

export Mesh

abstract AbstractMesh{T}
valuetype{T}(m::AbstractMesh{T}) = T

# UDim: dimension of the universe
# MDim: dimension of the mesh
# CType: type of the coordinates
# type Mesh{UDim,MDim,CType} <: AbstractMesh{CType}
#     vertices::Array{CType,2}
#     faces::Array{Int,2}
# end

# U: the dimension of the embedding space
# D1: the dimension of the mesh + 1 (!)
# T: the type of the coordinates of the vertices
type Mesh{U,D1,T} <: AbstractMesh{T}
    vertices::Array{Point{U,T},1}
    faces::Array{Vec{D1,Int},1}
end

vertextype{U,D1,T}(m::Mesh{U,D1,T}) = Point{U,T}

#vertices(m::Mesh) = m.vertices
#vertices(m::Mesh, V) = m.vertices[:,V]
function vertices(m::Mesh, I)
    r = zeros(vertextype(m), length(I))
    for i in I
        r[i] = m.vertices[i]
    end

    return r
end

#numvertices(m::Mesh) = size(m.vertices,2)
numvertices(m::Mesh) = length(m.vertices)
#numcells(m::Mesh) = size(m.faces,2)
numcells(m::Mesh) = length(m.faces)
#dimension(m::Mesh) = size(m.faces,1)-1
dimension{U,D1,T}(m::Mesh{U,D1,T}) = D1 - 1
#spacedimension(m::Mesh) = size(m.vertices,1)
spacedimension{U,D1,T}(m::Mesh{U,D1,T}) = U

function meshsegment{T<:Real}(length::T, delta::T, udim=2)

    num_segments = ceil(Int, length/delta)
    actual_delta = length / num_segments
    x = transpose(collect(0 : num_segments) * delta)
    y = zeros(x);

    vertices = zeros(T, udim, num_segments+1)
    vertices[1,:] = x
    vertices[2,:] = y

    faces = [ transpose(collect(1:num_segments));
              transpose(collect(2:num_segments+1))]

    return Mesh{udim,1,T}(vertices, faces)
end

function meshcircle{T<:Real}(radius::T, delta::T)

    circumf = 2 * pi *radius
    num_segments = ceil(Int, circumf / delta)
    delta = circumf / num_segments
    dalpha = delta / radius
    alpha = transpose(collect(0 : num_segments-1) * dalpha)
    vertices   = radius * [
        cos(alpha);
        sin(alpha) ]
    faces      = [
        transpose(1:num_segments);
        transpose(2:num_segments+1) ]
    faces[2,end] = 1

    return Mesh{2,1,T}(vertices, faces)
end

function meshrectangle{T}(width::T, height::T, delta::T, udim=3)

    @assert 2 <= udim

    UDim = 3
    MDim = 2

    nx = round(Int, ceil(width/delta));  nx = max(nx,1); dx = width/nx
    ny = round(Int, ceil(height/delta)); ny = max(ny,1); dy = height/ny

    xs = (0:nx) * dx
    ys = (0:ny) * dy

    #vertices = zeros(T, udim, (nx+1)*(ny+1))
    vertices = zeros(Point{UDim, T}, (nx+1)*(ny+1))
    k = 1
    for x in xs
        for y in ys
            p = zeros(T, UDim)
            p[1] = x
            p[2] = y
            vertices[k] = Point(p)
            k += 1
        end
    end

    #faces = zeros(Int, 3, 2*nx*ny)
    faces = zeros(Vec{3, Int}, 2*nx*ny)
    k = 1
    for i in 1 : nx
        for j in 1 : ny
            v11 = (i-1)*(ny+1) + j
            v12 = i*(ny+1) + j
            v21 = (i-1)*(ny+1) + (j+1)
            v22 = i*(ny+1) + (j+1)
            faces[k]   = Vec(v11,v21,v12)
            faces[k+1] = Vec(v21,v22,v12)
            k += 2
        end
    end

    #Mesh{UDim,MDim+1,T}(vertices, faces)
    Mesh(vertices, faces)
end


function meshfromfile(filename)
    open(filename) do f
        # multi-mesh files are not supported
        readline(f)

        # read the number of vertices and faces
        l = readline(f)

        sl = split(l)
        num_vertices = parse(Int, sl[1])
        num_faces    = parse(Int, sl[2])

        # read the vertex buffer
        vertices = zeros(Float64, 3, num_vertices)
        for i = 1 : num_vertices
            l = readline(f)
            vertices[:,i] = float(split(l))
        end

        # read the index buffer
        faces = zeros(Int64, 3, num_faces)
        for i = 1 : num_faces
            l = readline(f)
            #faces[:,i] = int(split(l))
            faces[:,i] = [parse(Int,s) for s in split(l)]
        end

        Mesh{3,2,Float64}(vertices, faces)

    end
end

# function to determine the boundary of a manifold
# Given a D-manifold, return a D-1 manifold representing
# the boundary and a array containing the indices of the
# vertices of the original mesh that corresond to he
# vertices of the boundary mesh.
function boundary{U,D1,T}(mesh::Mesh{U,D1,T})

    D = D1 - 1

    # vertices have no boundary
    @assert 0 < D

    # build a list of D-1 cells
    edges = cells(mesh, D-1)
    faces = cells(mesh, D)

    # get the edge-face connection matrix
    conn = connectivity(mesh, D-1, edges, faces)

    # find the edges that only have one edjacent face
    i = find(x -> x < 2, sum(abs(conn), 1))

    # create a mesh out of these
    bnd = Mesh{U,D1-1,T}(mesh.vertices, edges[i])
end

"""
This function takes an array of indices into the mesh vertex buffer
representing cells of a fixed dimension and returns an array of size
numvertices(mesh) x nmax integers containing in row i the indices
of the cells containing the i-th vertex. In addition the function returns
an array of length numvertices(mesh) containing the number of cells each
vertex is a member of.
"""
function vertextocellmap(mesh::Mesh, cells=mesh.faces)

    # first pass: determine the maximum number of cells
    # containing a given vertex

    numverts = numvertices(mesh)
    numcells = length(cells)
    numneighbors = zeros(Int, numverts)
    for i = 1 : numcells
        numneighbors[ cells[i] ] += 1
    end

    npos = -1
    vertstocells = fill(npos, numverts, maximum(numneighbors))
    numneighbors = zeros(Int, numverts)
    for i = 1 : numcells
        cell = cells[i]
        for j = 1 : length(cell)
            #v = cells[j,i]
            v = cell[j]
            k = (numneighbors[v] += 1)
            vertstocells[v,k] = i
        end
    end

    vertstocells, numneighbors
end

"""
Returns the cells of dimension k as an integer array of size k x N where N
is the number of cells of dimension k. The integers in this array represent
indices into the vertex buffer of the input mesh. Note that for the special
case k == 0 this function does not return any floating vertices.
"""
function cells(mesh::AbstractMesh, dim::Integer)

    meshdim = dimension(mesh)
    @assert 0 <= dim <= meshdim

    if dim == 0
        return [Vec(i) for i in 1:numvertices(mesh)]
    end

    if dim == meshdim return mesh.faces end

    C = numcells(mesh)
    simplices = zeros(Vec{dim+1,Int}, C*binomial(meshdim+1,dim+1))

    n = 1
    for c = 1 : C

        cell = mesh.faces[c]
        #cell = getsimplex(mesh, c)
        for simplex in combinations(cell,dim+1)

            simplices[n] = sort(simplex)
            n += 1
        end
    end

    simplices = unique(simplices)
end


function cells(f, mesh::AbstractMesh, k::Integer)

    @assert 0 <= k <= dimension(mesh)

    # compute upper bound for the number cells and allocate memory
    nc = numcells(mesh)
    kcells = zeros(Int, k+1, nc*binomial(dimension(mesh)+1, k+1))

    n = 1
    for c in 1 : nc

        dcell = mesh.faces[:,c]
        for kcell in combinations(dcell, k+1)

            if f(kcell)
                kcells[:,n] = sort(kcell)
                n += 1
            end
        end
    end

    kcells = kcells[:,1:n-1]
    kcells = unique(kcells,2)
end

function findfirst{N,T}(A::Array{T}, V::Union{Array{T,1},Vec{N,T}})
    I = zeros(Int, length(V))
    for (k,v) in enumerate(V)
        I[k] = Base.findfirst(A, v)
    end
    return I
end

# findfirst(Array, Vec)



# not exported
const levicivita_lut = cat(3,
    [0 0  0;  0 0 1; 0 -1 0],
    [0 0 -1;  0 0 0; 1  0 0],
    [0 1  0; -1 0 0; 0  0 0])

# Levi-Civita symbol of a permutation.
# The parity is computed by using the fact that a permutation is odd if and
# only if the number of even-length cycles is odd.
# Returns 1 is the permutarion is even, -1 if it is odd and 0 otherwise.
function levicivita{T<:Integer}(p::AbstractVector{T})
    n = length(p)

    if n == 3
        @inbounds valid = (0 < p[1] <= 3) * (0 < p[2] <= 3) * (0 < p[3] <= 3)
        return valid ? levicivita_lut[p[1], p[2], p[3]] : 0
    end

    todo = trues(n)
    first = 1
    cycles = flips = 0

    while cycles + flips < n
        first = findnext(todo, first)
        (todo[first] $= true) && return 0
        j = p[first]
        (0 < j <= n) || return 0
        cycles += 1
        while j ≠ first
            (todo[j] $= true) && return 0
            j = p[j]
            (0 < j <= n) || return 0
            flips += 1
        end
    end

    return iseven(flips) ? 1 : -1
end

# given a simplex and a face returns:
# +v if face is the v-th face of the simplex oriented according to the simplex
# -v if face is the v-th face of the simplex oriented oppositely to the simplex
# 0 is face is not a face of the simplex
function relorientation(face, simplex)

    v = setdiff(simplex, face)
    if length(v) != 1
        return 0
    end

    # find the position of the missing vertex
    v = v[1]
    i = Base.findfirst(simplex, v)
    s = (-1)^(i-1)

    # remove that vertex from the simplex
    face2 = Array(Int, length(simplex)-1)
    for j in 1 : i-1
        face2[j] = simplex[j]
    end
    for j in i : length(simplex)-1
        face2[j] = simplex[j+1]
    end
    #face2 = [ simplex[1:i-1]; simplex[i+1:end]]

    # get the permutation that maps face to face2
    p = findfirst(face2,face)

    return s * levicivita(p) * i
end

function connectivity(mesh::Mesh, k, kcells, mcells)

    if !( -1 <= k <= dimension(mesh))
        throw(ErrorException("k needs to be between zero and dim(mesh)"))
    end

    if k == -1
        return zeros(Int, numvertices(mesh), 0)
    end

    if k == dimension(mesh)
        return zeros(Int, 0, numcells(mesh))
    end

    @show size(kcells), typeof(kcells)
    @show size(mcells), typeof(mcells)

    vtok, _ = vertextocellmap(mesh, kcells)
    vtom, _ = vertextocellmap(mesh, mcells)

    const npos = -1

    #dimk = size(kcells,2)
    #dimm = size(mcells,2)
    dimk = length(kcells)
    dimm = length(mcells)

    @show dimk, dimm
    D = spzeros(Int, dimm, dimk)

    for v in 1 : numvertices(mesh)
        for i in vtok[v,:]
            i == npos && break
            kcell = kcells[i]
            for j in vtom[v,:]
                j == npos && break
                mcell = mcells[j]
                D[j,i] = sign(relorientation(kcell, mcell))
            end
        end
    end

    return D

end

# create pairs of volume cells from hypercells
# and the assumption of global orientation
# Every column of the returned array is a pair of indices
# If both ids are positive then they both refer to faces
# in the mesh structure. If the second number `id` is negative, this
# means that the coresponding edge lies on the boundary and that
# this edge has local index `-id` in the cell refererred to by the first index
function buildvolmpairs(mesh::Mesh, edges)

    @assert size(edges,1) == 2
    @assert ndims(edges) == 2

    numedges = size(edges, 2)
    facepairs = zeros(Int, 2, numedges)

    v2e, nn = vertextocellmap(mesh, mesh.faces)
    for e in 1:numedges

        edge = edges[:,e]

        # neighborhood of startvertex
        v = edge[1]
        n = nn[v]
        nbd1 = v2e[v,1:n]

        # neighborhood of endvertex
        v = edge[2]
        n = nn[v]
        nbd2 = v2e[v,1:n]

        # get the intersection
        nbd = intersect(nbd1, nbd2)

        n = length(nbd)
        @assert 0 <= n <= 2
        if n == 1 # boundary edge
            c = nbd[1]
            cell = mesh.faces[:,c]
            s = relorientation(edge, cell)
            facepairs[1,e] = c
            facepairs[2,e] = -abs(s)
        else # internal edge
            k = 0
            for c in nbd
                cell = mesh.faces[:,c]
                s = relorientation(edge, cell)
                if s > 0
                    facepairs[1,e] = c
                else
                    facepairs[2,e] = c
                end
            end
        end
    end

    facepairs

end