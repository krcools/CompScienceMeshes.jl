function meshsegment{T<:Real}(L::T, delta::T, udim=2)

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

  faces = Array{CT}(num_segments)
  for i in 1 : length(faces)
    faces[i] = SVector(i, i+1)
  end

  Mesh(vertices, faces)
end

function meshcircle{T<:Real}(radius::T, delta::T, udim=2)

  PT = SVector{udim,T}
  CT = SVector{2,Int}

    circumf = 2 * pi *radius
    num_segments = ceil(Int, circumf / delta)
    delta = circumf / num_segments
    dalpha = delta / radius
    alpha = collect(0 : num_segments-1) * dalpha

    vertices = Array{PT}(num_segments)
    for i in 1 : num_segments
        a = zeros(T, udim)
        a[1] = radius * cos(alpha[i])
        a[2] = radius * sin(alpha[i])
        vertices[i] = PT(a)
    end

    faces = Array{CT}(num_segments)
    for i in 1 : length(faces)-1
        faces[i] = SVector{2,Int}(i, i+1)
    end
    faces[end] = SVector{2,Int}(num_segments, 1)

    return Mesh(vertices, faces)
end


"""
    meshrectangle(width, height, delta, udim)

Create a mesh for a rectangle of width (along the x-axis) `width` and height (along
    the y-axis) `height`.

The target edge size is `delta` and the dimension of the
    embedding universe is `udim` (>= 2).
"""
function meshrectangle{T}(width::T, height::T, delta::T, udim=3)

  PT = SVector{udim,T}
  CT = SVector{3,Int}

    @assert 2 <= udim

    nx = round(Int, ceil(width/delta));  nx = max(nx,1); dx = width/nx
    ny = round(Int, ceil(height/delta)); ny = max(ny,1); dy = height/ny

    xs = (0:nx) * dx
    ys = (0:ny) * dy

    vertices = zeros(PT, (nx+1)*(ny+1))
    k = 1
    for x in xs
        for y in ys
            p = zeros(T, udim)
            p[1] = x
            p[2] = y
            vertices[k] = PT(p)
            k += 1
        end
    end

    faces = zeros(CT, 2*nx*ny)
    k = 1
    for i in 1 : nx
        for j in 1 : ny
            v11 = (i-1)*(ny+1) + j
            v12 = i*(ny+1) + j
            v21 = (i-1)*(ny+1) + (j+1)
            v22 = i*(ny+1) + (j+1)
            faces[k]   = CT(v11,v21,v12)
            faces[k+1] = CT(v21,v22,v12)
            k += 2
        end
    end

    Mesh(vertices, faces)
end


"""
    meshsphere(radius, delta)

Create a mesh of a sphere of radius `radius` by parsing a .geo script
    incorporating these parameters into the GMSH mesher.

The target edge size is `delta`.
"""
function meshsphere(radius, delta)
    s =
"""
lc = $delta;

Point(1)={0,0,0,lc};
Point(2)={$radius,0,0,lc};
Point(3)={0,$radius,0,lc};
Point(4)={-$radius,0,0,lc};
Point(5)={0,-$radius,0,lc};
Point(6)={0,0,$radius,lc};
Point(7)={0,0,-$radius,lc};

Circle(1)={2,1,3};
Circle(2)={3,1,4};
Circle(3)={4,1,5};
Circle(4)={5,1,2};
Circle(5)={6,1,3};
Circle(6)={3,1,7};
Circle(7)={7,1,5};
Circle(8)={5,1,6};

Line Loop(1)={-1,-6,-7,-4};
Line Loop(2)={1,-5,-8,4};
Line Loop(3)={-2,-3,7,6};
Line Loop(4)={2,3,8,5};

Ruled Surface(1)={1} In Sphere{1};
Ruled Surface(2)={2} In Sphere{1};
Ruled Surface(3)={3} In Sphere{1};
Ruled Surface(4)={4} In Sphere{1};
"""

    fn = tempname()
    io = open(fn, "w")
    try
        print(io, s)
    finally
        close(io)
    end

    # feed the file to gmsh
    fno = tempname()
    run(`gmsh $fn -2 -format msh -o $fno`)
    fdo = open(fno,"r")
    m = read_gmsh_mesh(fdo)
    close(fdo)
    rm(fno)
    rm(fn)

    return m

end


"""
    meshcuboid(width, height, length, delta)

Creates a mesh for a cuboid of width (along the x-axis) `width`, height (along
    the y-axis) `height` and length (along the z-axis) `length` by parsing a .geo script
    incorporating these parameters into the GMSH mesher.

The target edge size is `delta`.

"""
function meshcuboid(width, height, length, delta)
    s =
"""
lc = $delta;

Point(1)={0,0,0,lc};
Point(2)={0,0,$length,lc};
Point(3)={$width,0,$length,lc};
Point(4)={$width,0,0,lc};
Point(5)={0,$height,0,lc};
Point(6)={0,$height,$length,lc};
Point(7)={$width,$height,$length,lc};
Point(8)={$width,$height,0,lc};

Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};
Line(5)={5,6};
Line(6)={6,7};
Line(7)={7,8};
Line(8)={8,5};
Line(9)={1,5};
Line(10)={2,6};
Line(11)={3,7};
Line(12)={4,8};

Line Loop(1)={1,2,3,4};
Line Loop(2)={-1,9,5,-10};
Line Loop(3)={-2,10,6,-11};
Line Loop(4)={-3,11,7,-12};
Line Loop(5)={-4,12,8,-9};
Line Loop(6)={5,6,7,8};

Plane Surface(1)={1};
Plane Surface(2)={2};
Plane Surface(3)={3};
Plane Surface(4)={4};
Plane Surface(5)={5};
Plane Surface(6)={6};

Surface Loop(1)={1,2,3,4,5,6};

Volume(1)={1};

"""

    fn = tempname()
    io = open(fn, "w")
    try
        print(io, s)
    finally
        close(io)
    end

    # feed the file to gmsh
    fno = tempname()
    run(`gmsh $fn -2 -format msh -o $fno`)
    fdo = open(fno,"r")
    m = read_gmsh_mesh(fdo)
    close(fdo)
    rm(fno)
    rm(fn)

    return m

end


"""
    meshwaveguidepost(width, height, length, post_radius, delta)

Creates a mesh for a waveguide of width (along the x-axis) `width`, height (along
    the y-axis) `height` and length (along the z-axis) `length` by parsing a .geo script
    incorporating these parameters into the GMSH mesher.

A cylindrical post of radius `post_radius` is positioned vertically (i.e.
    parallel to y-axis) midway along the waveguide.

The target edge size is `delta`.
"""

function meshwaveguidepost(width, height, length, post_radius, delta)
    s =
"""

lc = $delta;

//--------------WAVEGUIDE POINTS--------------

Point(1) = {0.0, 0.0, 0.0, lc};
Point(2) = {0.0, $height, 0.0, lc};
Point(3) = {$width, $height, 0.0, lc};
Point(4) = {$width, 0.0, 0.0, lc};

Point(5) = {0.0, 0.0, $length, lc};
Point(6) = {0.0, $height, $length, lc};
Point(7) = {$width, $height, $length, lc};
Point(8) = {$width, 0.0, $length, lc};


//--------------POST POINTS-------------------


Point(9) = {($width)/2,0,($length)/2,lc};

Point(10) = {($width)/2,0,(($length)/2)-(($post_radius)/2),lc};
Point(11) = {(($width)/2)-(($post_radius)/2),0,($length)/2,lc};
Point(12) = {($width)/2,0,(($length)/2)+(($post_radius)/2),lc};
Point(13) = {(($width)/2)+(($post_radius)/2),0,($length)/2,lc};

Point(14) = {($width)/2,$height,($length)/2,lc};

Point(15) = {($width)/2,$height,(($length)/2)-(($post_radius)/2),lc};
Point(16) = {(($width)/2)-(($post_radius)/2),$height,($length)/2,lc};
Point(17) = {($width)/2,$height,(($length)/2)+(($post_radius)/2),lc};
Point(18) = {(($width)/2)+(($post_radius)/2),$height,($length)/2,lc};

//--------------WAVEGUIDE-POST CONTINUOUS SURFACE DEFINITION---------

//Part 1 - waveguide boundaries

Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};

Line(5)={5,6};
Line(6)={6,7};
Line(7)={7,8};
Line(8)={8,5};

Line(9)={1,5};
Line(10)={2,6};
Line(11)={3,7};
Line(12)={4,8};

//Part 2 - discontinuity boundaries

Circle(13) = {10,9,11};
Circle(14) = {11,9,12};
Circle(15) = {12,9,13};
Circle(16) = {13,9,10};

Circle(17) = {15,14,16};
Circle(18) = {16,14,17};
Circle(19) = {17,14,18};
Circle(20) = {18,14,15};

Line(21) = {10,15};
Line(22) = {11,16};
Line(23) = {12,17};
Line(24) = {13,18};

//Part 3 - creating skeleton

Line Loop(1)={-1,9,5,-10};

Line Loop(2)={-2,10,6,-11};
Line Loop(3)={17,18,19,20};

Line Loop(4)={-3,11,7,-12};

Line Loop(5)={-4,12,8,-9};
Line Loop(6)={-16,-15,-14,-13};

Line Loop(7) = {13,22,-17,-21};
Line Loop(8) = {14,23,-18,-22};
Line Loop(9) = {15,24,-19,-23};
Line Loop(10) = {16,21,-20,-24};

//Part 4 - creating surfaces

Plane Surface(1)={1};
Plane Surface(2)={3,2};
Plane Surface(3)={4};
Plane Surface(4)={6,5};

Ruled Surface(5) = {7};
Ruled Surface(6) = {8};
Ruled Surface(7) = {9};
Ruled Surface(8) = {10};

"""

    fn = tempname()
    io = open(fn, "w")
    try
        print(io, s)
    finally
        close(io)
    end

    # feed the file to gmsh
    fno = tempname()
    run(`gmsh $fn -2 -format msh -o $fno`)
    fdo = open(fno,"r")
    m = read_gmsh_mesh(fdo)
    # close(fdo)
    # rm(fno)
    # rm(fn)

    return m

end
