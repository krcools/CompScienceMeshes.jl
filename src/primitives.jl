using GmshTools


function meshgeo(geofile; physical=nothing, dim=2, tempname=tempname(), kwargs...)
    
    io = open(geofile)
    str = read(io, String)
    close(io)
    
    for (key,val) in kwargs
        key_str = String(key)
        pat = Regex("$key_str\\s*=\\s*[0-9]*.?[0-9]*;")
        sub = SubstitutionString("$key_str = $val;")
        str = replace(str, pat => sub)
    end

    println(str)
    # return
    
    temp_geo = tempname
    open(temp_geo, "w") do io
        print(io, str)
    end

    temp_msh = tempname * ".msh"
    gmsh.initialize()
    gmsh.option.setNumber("Mesh.MshFileVersion",2)
    gmsh.open(temp_geo)
    gmsh.model.mesh.generate(dim)
    gmsh.write(temp_msh)
    gmsh.finalize()

    # @show temp_msh

    # run(`gmsh $temp_geo -2 -format msh2 -o $temp_msh`)
    if dim == 2
        m = read_gmsh_mesh(temp_msh, physical=physical)
    elseif dim == 3
        m = read_gmsh3d_mesh(temp_msh, physical=physical)
    else
        error("gmsh files of dimension $(dim) cannot be read.")
    end
    
    rm(temp_msh)
    rm(temp_geo)
    
    return m
end

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


"""
    meshrectangle(width, height, delta, udim)

Create a mesh for a rectangle of width (along the x-axis) `width` and height (along
    the y-axis) `height`.

The target edge size is `delta` and the dimension of the
    embedding universe is `udim` (>= 2).

The mesh is oriented such that the normal is pointing down. This is subject to change.
"""
function meshrectangle(width::T, height::T, delta::T, udim=3; structured=true) where T
  if !structured
	  @assert udim==3 "Only 3D Unstructured mesh currently supported"
	  return meshrectangle_unstructured(width, height, delta)
  end

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
    meshsphere(;radius, h)

Create a mesh of a sphere of radius `radius` by parsing a .geo script
    incorporating these parameters into the GMSH mesher.

The target edge size is `delta`.
"""
function meshsphere(radius, delta; tempname=tempname())
    s = """
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

    fn = tempname
    io = open(fn, "w")
    try
        print(io, s)
    finally
        close(io)
    end
    fno = tempname * ".msh"

    gmsh.initialize()
    gmsh.option.setNumber("Mesh.MshFileVersion",2)
    gmsh.open(fn)
    gmsh.model.mesh.generate(2)
    gmsh.write(fno)
    gmsh.finalize()

    m = read_gmsh_mesh(fno)

    rm(fno)
    rm(fn)

    return m

end

meshsphere(;radius, h) = meshsphere(radius, h)

function meshsphere2(radius, h)
    fno = tempname() * ".msh"

    gmsh.initialize()
    gmsh.model.add("sphere")

    gmsh.model.geo.addPoint(0.0, 0.0, 0.0, h, 1)
    gmsh.model.geo.addPoint(radius, 0.0, 0.0, h, 2)
    gmsh.model.geo.addPoint(0.0, radius, 0.0, h, 3)
    gmsh.model.geo.addCircleArc(2, 1, 3, 1)
    gmsh.model.geo.addPoint(-radius, 0.0, 0.0, h, 4)
    gmsh.model.geo.addPoint(0.0, -radius, 0.0, h, 5)
    gmsh.model.geo.addCircleArc(3, 1, 4, 2)
    gmsh.model.geo.addCircleArc(4, 1, 5, 3)
    gmsh.model.geo.addCircleArc(5, 1, 2, 4)
    gmsh.model.geo.addPoint(0.0, 0.0, -radius, h, 6)
    gmsh.model.geo.addPoint(0.0, 0.0, radius, h, 7)
    gmsh.model.geo.addCircleArc(3, 1, 6, 5)
    gmsh.model.geo.addCircleArc(6, 1, 5, 6)
    gmsh.model.geo.addCircleArc(5, 1, 7, 7)
    gmsh.model.geo.addCircleArc(7, 1, 3, 8)
    gmsh.model.geo.addCircleArc(2, 1, 7, 9)
    gmsh.model.geo.addCircleArc(7, 1, 4, 10)
    gmsh.model.geo.addCircleArc(4, 1, 6, 11)
    gmsh.model.geo.addCircleArc(6, 1, 2, 12)
    gmsh.model.geo.addCurveLoop([2, 8, -10], 13)
    gmsh.model.geo.addSurfaceFilling([13], 14)
    gmsh.model.geo.addCurveLoop([10, 3, 7], 15)
    gmsh.model.geo.addSurfaceFilling([15], 16)
    gmsh.model.geo.addCurveLoop([-8, -9, 1], 17)
    gmsh.model.geo.addSurfaceFilling([17], 18)
    gmsh.model.geo.addCurveLoop([-11, -2, 5], 19)
    gmsh.model.geo.addSurfaceFilling([19], 20)
    gmsh.model.geo.addCurveLoop([-5, -12, -1], 21)
    gmsh.model.geo.addSurfaceFilling([21], 22)
    gmsh.model.geo.addCurveLoop([-3, 11, 6], 23)
    gmsh.model.geo.addSurfaceFilling([23], 24)
    gmsh.model.geo.addCurveLoop([-7, 4, 9], 25)
    gmsh.model.geo.addSurfaceFilling([25], 26)
    gmsh.model.geo.addCurveLoop([-4, 12, -6], 27)
    gmsh.model.geo.addSurfaceFilling([27], 28)
    gmsh.model.geo.addSurfaceLoop([28, 26, 16, 14, 20, 24, 22, 18], 29)
    gmsh.model.geo.addVolume([29], 30)
    gmsh.model.geo.addPhysicalGroup(2, [28, 26, 16, 14, 20, 24, 22, 18], 1)
    gmsh.model.geo.addPhysicalGroup(3, [30], 2)

    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.MshFileVersion",2)
    gmsh.model.mesh.generate(2)
    gmsh.fltk.run()
    gmsh.write(fno)
    gmsh.finalize()

    m = read_gmsh_mesh(fno)
    rm(fno)
    return m

end

meshsphere2(;radius, h) = meshsphere2(radius, h)

"""
    meshtorus(innerradius, outerradius, h)
    meshtorus(;innerradius, outerradius, h)

Create a mesh of a torus of 2 radii `innerradius` and `outerradius`

The target edge size is `h`.
"""
function meshtorus(innerradius, outerradius, h)
    @assert innerradius < outerradius
    
    fno = tempname() * ".msh"
    center = (outerradius + innerradius)/2
    radius = (outerradius - innerradius)/2

    gmsh.initialize()
    gmsh.model.add("torus")
    gmsh.model.geo.addPoint(0.0, 0.0, 0.0, h, 1)
    gmsh.model.geo.addPoint(0.0, 0.0, radius, h, 2)
    gmsh.model.geo.addPoint(0.0, 0.0, -radius, h, 3)

    # 1st quarter
    gmsh.model.geo.addPoint(center, 0.0, 0.0, h, 4)
    gmsh.model.geo.addPoint(outerradius, 0.0, 0.0, h, 5)
    gmsh.model.geo.addPoint(center, 0.0, radius, h, 6)
    gmsh.model.geo.addPoint(innerradius, 0.0, 0.0, h, 7)
    gmsh.model.geo.addPoint(center, 0.0, -radius, h, 8)
    gmsh.model.geo.addCircleArc(5, 4, 6, 1)
    gmsh.model.geo.addCircleArc(6, 4, 7, 2)
    gmsh.model.geo.addCircleArc(7, 4, 8, 3)
    gmsh.model.geo.addCircleArc(8, 4, 5, 4)

    gmsh.model.geo.addPoint(0.0, center, 0.0, h, 9)
    gmsh.model.geo.addPoint(0.0, outerradius, 0.0, h, 10)
    gmsh.model.geo.addPoint(0.0, center, radius, h, 11)
    gmsh.model.geo.addPoint(0.0, innerradius, 0.0, h, 12)
    gmsh.model.geo.addPoint(0.0, center, -radius, h, 13)
    gmsh.model.geo.addCircleArc(10, 9, 11, 5)
    gmsh.model.geo.addCircleArc(11, 9, 12, 6)
    gmsh.model.geo.addCircleArc(12, 9, 13, 7)
    gmsh.model.geo.addCircleArc(13, 9, 10, 8)
    
    gmsh.model.geo.addCircleArc(5, 1, 10, 9)
    gmsh.model.geo.addCircleArc(6, 2, 11, 10)
    gmsh.model.geo.addCircleArc(7, 1, 12, 11)
    gmsh.model.geo.addCircleArc(8, 3, 13, 12)

    gmsh.model.geo.addCurveLoop([1, 10, -5, -9], 13)
    gmsh.model.geo.addSurfaceFilling([13], 1)
    gmsh.model.geo.addCurveLoop([2, 11, -6, -10], 14)
    gmsh.model.geo.addSurfaceFilling([14], 2)
    gmsh.model.geo.addCurveLoop([3, 12, -7, -11], 15)
    gmsh.model.geo.addSurfaceFilling([15], 3)
    gmsh.model.geo.addCurveLoop([4, 9, -8, -12], 16)
    gmsh.model.geo.addSurfaceFilling([16], 4)

    sf = gmsh.model.geo.copy([(2, 1), (2, 2), (2, 3), (2, 4)])
    gmsh.model.geo.rotate(sf, 0, 0, 0, 0, 0, 1, -pi/2)
    sf = gmsh.model.geo.copy([(2, 1), (2, 2), (2, 3), (2, 4)])
    gmsh.model.geo.rotate(sf, 0, 0, 0, 0, 0, 1, pi/2)
    sf = gmsh.model.geo.copy([(2, 1), (2, 2), (2, 3), (2, 4)])
    gmsh.model.geo.rotate(sf, 0, 0, 0, 0, 0, 1, pi)

    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.MshFileVersion",2)
    gmsh.model.mesh.generate(2)
    # gmsh.fltk.run()
    gmsh.write(fno)
    gmsh.finalize()

    m = read_gmsh_mesh(fno)
    rm(fno)
    return m
end

meshtorus(;innerradius, outerradius, h) = meshtorus(innerradius, outerradius, h)

"""
not working yet
"""
function tetmeshsphere(radius,delta; tempname=tempname())
    s = """
    lc = $delta;

    Point(1)={0,0,0,lc};
    Point(2)={$radius,0,0,lc};
    Point(3)={-$radius,0,0,lc};
    Point(4)={0,$radius,0,lc};
    Point(5)={0,-$radius,0,lc};
    Point(6)={0,0,$radius,lc};
    Point(7)={0,0,-$radius,lc};

    Circle(1)={7,1,3};
    Circle(2)={3,1,6};
    Circle(3)={6,1,2};
    Circle(4)={2,1,7};
    Circle(5)={7,1,4};
    Circle(6)={4,1,6};
    Circle(7)={6,1,5};
    Circle(8)={5,1,7};
    Circle(9)={3,1,5};
    Circle(10)={5,1,2};
    Circle(11)={2,1,4};
    Circle(12)={4,1,3};

    Curve Loop(1)={8,-4,-10};
    Curve Loop(2)={11,-5,-4};
    Curve Loop(3)={5,12,-1};
    Curve Loop(4)={1,9,8};
    Curve Loop(5)={6,-2,-12};
    Curve Loop(6)={9,-7,-2};
    Curve Loop(7)={7,10,-3};
    Curve Loop(8)={6,3,11};

    Surface(1)={1};
    Surface(2)={2};
    Surface(3)={3};
    Surface(4)={4};
    Surface(5)={5};
    Surface(6)={6};
    Surface(7)={7};
    Surface(8)={8};

    Surface Loop(1) = {4, 3, 2, 8, 5, 6, 7, 1};
    Volume(1) = {1};
    Physical Volume(1) = {1};
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
    gmsh.option.setNumber("Mesh.MshFileVersion",2)
    gmsh.open(fn)
    gmsh.model.mesh.generate(3)
    gmsh.write(fno)
    gmsh.finalize()

    m = read_gmsh3d_mesh(fno)

    rm(fno)
    rm(fn)

    return m
end



function meshball(;radius, h, tempname=tempname())

    fn = joinpath(@__DIR__,"geos/structured_sphere.geo")
    io = open(fn)
    str = read(io, String)
    close(io)

    str = replace(str, "lc = 1.0;" => "lc = $h;")
    str = replace(str, "r = 1.0;" => "r = $radius;")

    temp_geo = tempname
    open(temp_geo, "w") do io
        print(io, str)
    end

    temp_msh = tempname * ".msh"
    gmsh.initialize()
    gmsh.option.setNumber("Mesh.MshFileVersion",2)
    gmsh.open(temp_geo)
    gmsh.model.mesh.generate(3)
    gmsh.write(temp_msh)
    gmsh.finalize()

    m = read_gmsh3d_mesh(temp_msh)

    rm(temp_msh)
    rm(temp_geo)

    return m
end


function meshcylinder(;radius, height, h, tempname=tempname())

    fn = joinpath(@__DIR__,"geos/cylinder3.geo")
    io = open(fn)
    str = read(io, String)
    close(io)

    str = replace(str, "r = 1.0;" => "r = $radius;")
    str = replace(str, "z = 1.0;" => "z = $height;")
    str = replace(str, "h = 1.0;" => "h = $h;")

    temp_geo = tempname
    open(temp_geo, "w") do io
        print(io, str)
    end

    temp_msh = tempname * ".msh"
    gmsh.initialize()
    gmsh.option.setNumber("Mesh.MshFileVersion",2)
    gmsh.open(temp_geo)
    gmsh.model.mesh.generate(3)
    gmsh.write(temp_msh)
    gmsh.finalize()
    m = read_gmsh3d_mesh(temp_msh)

    rm(temp_msh)
    rm(temp_geo)

    return m
end


"""
    meshcuboid(width, height, length, delta)

Creates a mesh for a cuboid of width (along the x-axis) `width`, height (along
    the y-axis) `height` and length (along the z-axis) `length` by parsing a .geo script
    incorporating these parameters into the GMSH mesher.

The target edge size is `delta`.
physical => in {"TopPlate", "BottomPlate", "SidePlates", "OpenBox"} extracts and
returns only the specified part of the cuboid

"""
function meshcuboid(width::T, height::T, length::T, delta::T;physical="ClosedBox", tempname=tempname()) where T
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

Line Loop(1)={-1,-2,-3,-4};
Line Loop(2)={1,-9,-5,10};
Line Loop(3)={2,-10,-6,11};
Line Loop(4)={3,-11,-7,12};
Line Loop(5)={4,-12,-8,9};
Line Loop(6)={5,6,7,8};

Plane Surface(1)={1};
Plane Surface(2)={2};
Plane Surface(3)={3};
Plane Surface(4)={4};
Plane Surface(5)={5};
Plane Surface(6)={6};

//classify parts of geometry
Physical Surface("TopPlate") = {3};
Physical Surface("BottomPlate") = {5};
Physical Surface("SidePlates") = {1,2,4,6};
Physical Surface("OpenBox") = {1,2,4,5,6};
Physical Surface("ClosedBox") = {1,2,3,4,5,6};

Surface Loop(1)={1,2,3,4,5,6};
Volume(1)={1};

"""

    fn = tempname
    io = open(fn, "w")
    try
        print(io, s)
    finally
        close(io)
    end

    # feed the file to gmsh
    fno = tempname * ".msh"

    gmsh.initialize()
    gmsh.option.setNumber("Mesh.MshFileVersion",2)
    gmsh.open(fn)
    gmsh.model.mesh.generate(2)
    gmsh.write(fno)
    gmsh.finalize()

    # m = read_gmsh_mesh(fno)
    m = read_gmsh_mesh(fno,physical=physical,T=T)

    rm(fno)
    rm(fn)

    return m

end

"""
    meshsquaretorus(width, height, holewidth, h)
    meshsquaretorus(;width, height, holewidth, h)

Create a mesh of a square torus of size `width` and `height` with a hole of size `holewidth` and `height`

The target edge size is `h`.
"""
function meshsquaretorus(width, height, holewidth, h)
    @assert holewidth < width
    
    fno = tempname() * ".msh"
    gmsh.initialize()
    gmsh.model.add("squaretorus")

    # bottom plate
    gmsh.model.geo.addPoint(holewidth/2, -holewidth/2, -height/2, h, 1)
    gmsh.model.geo.addPoint(holewidth/2, holewidth/2, -height/2, h, 2)
    gmsh.model.geo.addPoint(-holewidth/2, holewidth/2, -height/2, h, 3)
    gmsh.model.geo.addPoint(-holewidth/2, -holewidth/2, -height/2, h, 4)
    gmsh.model.geo.addLine(1, 2, 1)
    gmsh.model.geo.addLine(2, 3, 2)
    gmsh.model.geo.addLine(3, 4, 3)
    gmsh.model.geo.addLine(4, 1, 4)
    gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 101)
    
    gmsh.model.geo.addPoint(width/2, -width/2, -height/2, h, 5)
    gmsh.model.geo.addPoint(width/2, width/2, -height/2, h, 6)
    gmsh.model.geo.addPoint(-width/2, width/2, -height/2, h, 7)
    gmsh.model.geo.addPoint(-width/2, -width/2, -height/2, h, 8)
    gmsh.model.geo.addLine(5, 6, 5)
    gmsh.model.geo.addLine(6, 7, 6)
    gmsh.model.geo.addLine(7, 8, 7)
    gmsh.model.geo.addLine(8, 5, 8)

    gmsh.model.geo.addCurveLoop([5, 6, 7, 8], 102)
    gmsh.model.geo.addPlaneSurface([101, 102], 1)
    gmsh.model.geo.extrude([(2, 1)], 0, 0, height)

    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.MshFileVersion",2)
    gmsh.model.mesh.generate(2)
    # gmsh.fltk.run()
    gmsh.write(fno)
    gmsh.finalize()

    m = read_gmsh_mesh(fno)
    rm(fno)
    return m
end

meshsquaretorus(;width, height, holewidth, h) = meshsquaretorus(width, height, holewidth, h)

"""
    meshsquaretorus4(width, height, holewidth, h)
    meshsquaretorus4(;width, height, holewidth, h)

Create a mesh of a square torus of size `width` and `height` with 4 holes of size `holewidth` and `height`

The target edge size is `h`.
"""
function meshsquaretorus4(width, height, holewidth, h)
    @assert (2*holewidth) < width
    
    fno = tempname() * ".msh"
    gmsh.initialize()
    gmsh.model.add("squaretorus4")

    # bottom plate
    gmsh.model.geo.addPoint(width/2, -width/2, -height/2, h, 1)
    gmsh.model.geo.addPoint(width/2, width/2, -height/2, h, 2)
    gmsh.model.geo.addPoint(-width/2, width/2, -height/2, h, 3)
    gmsh.model.geo.addPoint(-width/2, -width/2, -height/2, h, 4)
    gmsh.model.geo.addLine(1, 2, 1)
    gmsh.model.geo.addLine(2, 3, 2)
    gmsh.model.geo.addLine(3, 4, 3)
    gmsh.model.geo.addLine(4, 1, 4)
    gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 101)

    gmsh.model.geo.addPoint(holewidth/2 + width/4, -holewidth/2 - width/4, -height/2, h, 5)
    gmsh.model.geo.addPoint(holewidth/2 + width/4, holewidth/2 - width/4, -height/2, h, 6)
    gmsh.model.geo.addPoint(-holewidth/2 + width/4, holewidth/2 - width/4, -height/2, h, 7)
    gmsh.model.geo.addPoint(-holewidth/2 + width/4, -holewidth/2 - width/4, -height/2, h, 8)
    gmsh.model.geo.addLine(5, 6, 5)
    gmsh.model.geo.addLine(6, 7, 6)
    gmsh.model.geo.addLine(7, 8, 7)
    gmsh.model.geo.addLine(8, 5, 8)

    gmsh.model.geo.addCurveLoop([5, 6, 7, 8], 102)

    hole2 = gmsh.model.geo.copy([(1, 5), (1, 6), (1, 7), (1, 8)])
    gmsh.model.geo.translate(hole2, 0, width/2, 0)
    curveloop2 = [tag for (dim, tag) in hole2]
    gmsh.model.geo.addCurveLoop(curveloop2, 103)

    hole3 = gmsh.model.geo.copy([(1, 5), (1, 6), (1, 7), (1, 8)])
    gmsh.model.geo.translate(hole3, -width/2, width/2, 0)
    curveloop3 = [tag for (dim, tag) in hole3]
    gmsh.model.geo.addCurveLoop(curveloop3, 104)

    hole4 = gmsh.model.geo.copy([(1, 5), (1, 6), (1, 7), (1, 8)])
    gmsh.model.geo.translate(hole4, -width/2, 0, 0)
    curveloop4 = [tag for (dim, tag) in hole4]
    gmsh.model.geo.addCurveLoop(curveloop4, 105)

    gmsh.model.geo.addPlaneSurface([101, 102, 103, 104, 105], 1)
    gmsh.model.geo.extrude([(2, 1)], 0, 0, height)

    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.MshFileVersion",2)
    gmsh.model.mesh.generate(2)
    # gmsh.fltk.run()
    gmsh.write(fno)
    gmsh.finalize()

    m = read_gmsh_mesh(fno)
    rm(fno)
    return m
end

meshsquaretorus4(;width, height, holewidth, h) = meshsquaretorus4(width, height, holewidth, h)

"""
    meshrectangle_unstructured(width, height, delta)
	Meshes unstructured rectangle (Delaunay Triangulation)
"""
function meshrectangle_unstructured(width, height, delta; tempname=tempname())
    s =
		"""
		lc = $delta;

		Point(1)={0,0,0,lc};
		Point(4)={$width,0,0,lc};
		Point(5)={0,$height,0,lc};
		Point(8)={$width,$height,0,lc};

		Line(4)={4,1};
		Line(8)={8,5};
		Line(9)={1,5};
		Line(12)={4,8};

		Line Loop(5)={4,-12,-8,9};

		Plane Surface(5)={5};
		"""

    fn = tempname
    io = open(fn, "w")
    try
        print(io, s)
    finally
        close(io)
    end

    # feed the file to gmsh
    fno = tempname * ".msh"

    gmsh.initialize()
    gmsh.option.setNumber("Mesh.MshFileVersion",2)
    gmsh.open(fn)
    gmsh.model.mesh.generate(2)
    gmsh.write(fno)
    gmsh.finalize()

    m = read_gmsh_mesh(fno)

    rm(fno)
    rm(fn)
    return m

end


function meshmobius(;h)

    m = meshrectangle(2pi, 2.0, h, 3)
    m = translate(m, point(-pi, -1, 0))
    for (i,v) in enumerate(m.vertices)
        s,t = v[1], v[2]
        x = 2*point(cos(s), sin(s), 0) + t*(cos(s/2)*point(cos(s), sin(s),0) + sin(s/2)*point(0,0,1))
        m.vertices[i] = x
    end

    return m
end