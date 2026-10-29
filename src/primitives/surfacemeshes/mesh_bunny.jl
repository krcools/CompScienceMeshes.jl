using Downloads

const bunny_url = "https://raw.githubusercontent.com/alecjacobson/common-3d-test-models/master/data/stanford-bunny.obj"

"""
    meshbunny(h; order=1, size=1.0, obj=nothing)

Create a closed mesh of the Stanford bunny with target edge size `h`. The bunny is
centred at the origin, scaled to a largest dimension of `size`, and its ears point
in the ``y``-direction. The mesh is oriented with outward pointing normals.

With `order=1` a flat triangular mesh is returned. With `order=3` the result is a
curvilinear mesh of cubic triangles: every flat triangle is replaced by a PN triangle
(Vlachos et al., Curved PN triangles, ACM I3D 2001) that interpolates the vertices
and the normals of the scanned surface there.

The geometry is the range scan from the Stanford 3D Scanning Repository
(https://graphics.stanford.edu/data/3Dscanrep/). Its five holes in the base are closed
before the scan is remeshed by gmsh. The scan is downloaded on first use and cached in
the Julia depot; pass the path of a local copy of `stanford-bunny.obj` as `obj` to
work offline.
"""
function meshbunny(h; order=1, size=1.0, obj=nothing)
    obj === nothing && (obj = bunny_scanfile())
    V, F = bunny_scan(obj, size)

    fno = tempname() * ".msh"
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    # Hand the triangulation to gmsh as a discrete surface (keeps the connectivity exact)
    used = sort(unique(reduce(vcat, F)))
    tag = Dict(v => i for (i, v) in enumerate(used))
    s = gmsh.model.addDiscreteEntity(2)
    gmsh.model.mesh.addNodes(2, s, collect(1:length(used)), reduce(vcat, Vector.(V[used])))
    gmsh.model.mesh.addElementsByType(s, 2, Int[], [tag[v] for f in F for v in f])
    # Split the triangulation into patches gmsh can reparametrise (cf. gmsh tutorial t13).
    # The bunny has no sharp features, so split on topology only: smaller angles cut the
    # noisy scan into patches with invalid boundaries.
    gmsh.model.mesh.classifySurfaces(π, true, true, π)
    gmsh.model.mesh.createGeometry()
    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMin", h)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2)
    gmsh.model.mesh.generate(2)
    gmsh.write(fno)
    gmsh.finalize()
    m = bunny_orient(read_gmsh_mesh(fno))
    rm(fno)

    order == 1 && return m
    order == 3 || throw(ArgumentError("meshbunny supports order 1 and 3 only."))
    # The high order nodes gmsh places itself (setOrder) follow the noise of the scan
    # and fold about one in ten cubic elements; PN triangles are smooth by construction.
    return bunny_pntriangles(m, bunny_normals(m, V, F, h / 2))
end

function bunny_scanfile()
    dir = joinpath(first(DEPOT_PATH), "scratchspaces", "3e66a162-7b8c-5da0-b8f8-124ecd2c3ae1", "bunny")
    fn = joinpath(dir, "stanford-bunny.obj")
    if !isfile(fn)
        mkpath(dir)
        @info "Downloading the Stanford bunny to $fn"
        Downloads.download(bunny_url, fn)
    end
    return fn
end

# Watertight scan, scaled to a largest dimension of `size` and centred at the origin
function bunny_scan(fn, size)
    V = SVector{3,Float64}[]
    F = Vector{Int}[]
    for l in eachline(fn)
        s = split(l)
        isempty(s) && continue
        s[1] == "v" && push!(V, SVector{3,Float64}(parse.(Float64, s[2:4])))
        s[1] == "f" && push!(F, parse.(Int, first.(split.(s[2:4], "/"))))
    end
    bunny_fillholes!(V, F)
    lo, hi = reduce((a, b) -> min.(a, b), V), reduce((a, b) -> max.(a, b), V)
    V = [(v - (lo + hi) / 2) * size / maximum(hi - lo) for v in V]
    return V, F
end

# Close every hole by a fan of triangles from the hole's centroid. Fan triangles
# traverse each boundary edge opposite to the triangle it belongs to.
function bunny_fillholes!(V, F)
    count = Dict{Tuple{Int,Int},Int}()
    for f in F, k in 1:3
        e = minmax(f[k], f[mod1(k+1, 3)])
        count[e] = get(count, e, 0) + 1
    end
    next = Dict{Int,Int}()
    for f in F, k in 1:3
        a, b = f[k], f[mod1(k+1, 3)]
        count[minmax(a, b)] == 1 && (next[b] = a)
    end
    while !isempty(next)
        v0 = first(keys(next))
        loop = [v0]
        while (v = next[loop[end]]) != v0
            push!(loop, v)
        end
        foreach(v -> delete!(next, v), loop)
        push!(V, sum(V[loop]) / length(loop))
        c = length(V)
        for k in eachindex(loop)
            push!(F, [loop[k], loop[mod1(k+1, length(loop))], c])
        end
    end
    return V, F
end

# Make the orientations of the triangles of a closed surface consistent by propagating
# across shared edges, then flip everything if the normals point inward.
function bunny_orient(m)
    faces = [collect(c.indices) for c in cells(m)]
    edge2faces = Dict{Tuple{Int,Int},Vector{Int}}()
    for (f, c) in enumerate(faces), k in 1:3
        a, b = c[k], c[mod1(k+1, 3)]
        push!(get!(edge2faces, minmax(a, b), Int[]), f)
    end
    @assert all(length(v) == 2 for v in values(edge2faces)) "mesh is not a closed manifold"
    hasedge(c, a, b) = any(c[k] == a && c[mod1(k+1, 3)] == b for k in 1:3)
    done = falses(length(faces))
    done[1] = true
    stack = [1]
    while !isempty(stack)
        f = pop!(stack)
        c = faces[f]
        for k in 1:3
            a, b = c[k], c[mod1(k+1, 3)]
            for g in edge2faces[minmax(a, b)]
                (g == f || done[g]) && continue
                # a consistently oriented neighbour traverses the shared edge as b → a
                hasedge(faces[g], a, b) && reverse!(faces[g])
                done[g] = true
                push!(stack, g)
            end
        end
    end
    @assert all(done) "mesh is not connected"
    v = vertices(m)
    vol = sum(dot(v[c[1]], cross(v[c[2]], v[c[3]])) for c in faces) / 6
    vol < 0 && foreach(reverse!, faces)
    return Mesh(v, [SimplexGraph(index(c...)) for c in faces])
end

# Unit normals at the vertices of `m`, averaged (area weighted) over the faces of
# the dense scan (V, F) whose centroids lie within distance `r`, signed outward like `m`.
function bunny_normals(m, V, F, r)
    ctr = [sum(V[f]) / 3 for f in F]
    nrm = [cross(V[f[2]] - V[f[1]], V[f[3]] - V[f[1]]) / 2 for f in F]
    bin(p) = Tuple(floor.(Int, p / r))
    grid = Dict{NTuple{3,Int},Vector{Int}}()
    for (i, c) in enumerate(ctr)
        push!(get!(grid, bin(c), Int[]), i)
    end
    vn = zeros(SVector{3,Float64}, numvertices(m))    # outward normals of the coarse mesh
    for c in cells(m)
        a, b, d = vertices(m)[c.indices]
        foreach(v -> vn[v] += cross(b - a, d - a), c.indices)
    end
    map(enumerate(vertices(m))) do (i, p)
        b = bin(p)
        n = sum((nrm[j] for k in Iterators.product(-1:1, -1:1, -1:1) for j in get(grid, b .+ k, Int[])
            if norm(ctr[j] - p) < r); init=zero(SVector{3,Float64}))
        n = normalize(n)
        dot(n, vn[i]) < 0 ? -n : n
    end
end

# Cubic Bézier PN triangle with corners P and corner normals N, as a function of the
# barycentric coordinates (l1, l2, l3).
function bunny_pnpatch(P, N)
    b(i, j) = (2P[i] + P[j] - dot(P[j] - P[i], N[i]) * N[i]) / 3    # control point on edge i→j next to i
    b210, b120, b021, b012, b102, b201 = b(1,2), b(2,1), b(2,3), b(3,2), b(3,1), b(1,3)
    E = (b210 + b120 + b021 + b012 + b102 + b201) / 6
    b111 = E + (E - sum(P) / 3) / 2
    return (l1, l2, l3) -> P[1]*l1^3 + P[2]*l2^3 + P[3]*l3^3 +
        3b210*l1^2*l2 + 3b120*l1*l2^2 + 3b021*l2^2*l3 + 3b012*l2*l3^2 +
        3b102*l1*l3^2 + 3b201*l1^2*l3 + 6b111*l1*l2*l3
end

# Cubic curvilinear mesh interpolating PN triangles through the cells of the flat mesh
# `m` with vertex normals `vn`. Nodes follow gmsh's 10-node triangle ordering; the
# geometry of an edge only depends on its end points, so neighbouring cells conform.
function bunny_pntriangles(m, vn)
    verts = copy(vertices(m))
    edgenode = Dict{Tuple{Int,Int},Int}()    # (a, b) => node on edge ab at 1/3 from a
    faces = map(cells(m)) do c
        v = c.indices
        f = bunny_pnpatch(verts[v], vn[v])
        node(i, j, l) = get!(edgenode, (v[i], v[j])) do
            push!(verts, f(l...))
            length(verts)
        end
        e = (node(1, 2, (2/3, 1/3, 0)), node(2, 1, (1/3, 2/3, 0)),
             node(2, 3, (0, 2/3, 1/3)), node(3, 2, (0, 1/3, 2/3)),
             node(3, 1, (1/3, 0, 2/3)), node(1, 3, (2/3, 0, 1/3)))
        push!(verts, f(1/3, 1/3, 1/3))
        HOSimplexGraph{3,3,10}(SVector(v[1], v[2], v[3], e..., length(verts)))
    end
    return CurvilinearMesh(verts, faces)
end
