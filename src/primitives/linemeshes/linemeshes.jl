"""
    meshsegment(L::T, delta::T, udim=2) where T<:Real

"""
function meshsegment(L::T, delta::T, udim=2) where T<:Real

    PT = SVector{udim, T}
    CT = SVector{2,Int}

    num_segments = ceil(Int, L/delta)
    actual_delta = L/num_segments
    x = collect(0:num_segments) * actual_delta

    vertices = zeros(PT, num_segments+1)
    for i in 1 : length(vertices)
        a = zeros(T, udim)
        a[1] = x[i]
        vertices[i] = PT(a)
    end

    faces = Array{CT}(undef,num_segments)
    for i in 1 : length(faces)
        faces[i] = SVector(i, i+1)
    end

    Mesh(vertices, faces)
end
  
"""
    meshcircle(radius::T, delta::T, udim=2) where T<:Real

"""
function meshcircle(radius::T, delta::T, udim=2) where T<:Real

PT = SVector{udim,T}
CT = SVector{2,Int}

    circumf = 2 * pi *radius
    num_segments = ceil(Int, circumf / delta)
    delta = circumf / num_segments
    dalpha = delta / radius
    alpha = collect(0 : num_segments-1) * dalpha

    vertices = Array{PT}(undef,num_segments)
    for i in 1 : num_segments
        a = zeros(T, udim)
        a[1] = radius * cos(alpha[i])
        a[2] = radius * sin(alpha[i])
        vertices[i] = PT(a)
    end

    faces = Array{CT}(undef,num_segments)
    for i in 1 : length(faces)-1
        faces[i] = SVector{2,Int}(i, i+1)
    end
    faces[end] = SVector{2,Int}(num_segments, 1)

    return Mesh(vertices, faces)
end

function gmshcircle(radius::T, delta::T, tempname=tempname(); order::Int=1) where T<:Real
    s = """
lc = $delta;

Point(1) = {0, 0, 0, lc};
Point(2) = {$radius, 0, 0, lc};
Point(3) = {0, $radius, 0, lc};
Point(4) = {-$radius, 0, 0, lc};
Point(5) = {0, -$radius, 0, lc};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

Physical Curve("boundary") = {1, 2, 3, 4};
"""

    fn = tempname
    io = open(fn, "w")
    try
        print(io, s)
    finally
        close(io)
    end

    fno = tempname * ".msh"

    gmsh.initialize()
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.0)
    gmsh.open(fn)

    # Choose element order (1: linear segments, 2: quadratic)
    gmsh.option.setNumber("Mesh.ElementOrder", order)

    # Generate only the 1D mesh (curves)
    gmsh.model.mesh.generate(1)

    # If you *don’t* define Physicals, enable the next line instead:
    # gmsh.option.setNumber("Mesh.SaveAll", 1)

    gmsh.write(fno)
    gmsh.finalize()

    m = load_gmsh_mesh(fno,
        element=:line,
        udim=2,
        vertextype=Float64,
        order=order,
        sort=false
    )

    #rm(fno)
    rm(fn)

    return m
end