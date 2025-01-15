using GmshTools

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

"""
    meshball(;radius, h, tempname=tempname())

"""
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

